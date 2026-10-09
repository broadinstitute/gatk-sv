"""Verify single-record manta/dragen tloc auto-resolution in svtk resolve.

Semantics under test (mirrors the two-mate path's contract):
  - CANDIDATE_SINGLE_TLOC resolution is OPT-IN (resolve_single_tlocs=True,
    set by the --resolve-single-tlocs CLI flag, passed only by
    mantatloccheck.sh); with the flag off a single manta BND stays an
    unresolved SINGLE_ENDER exactly as on main;
  - single interchromosomal manta/dragen BND with CHR2/END2 + PE evidence
    resolves to SVTYPE=CTX;
  - the CPX_TYPE label comes from the cytoband arms (same arms -> PP/QQ,
    different arms -> PQ/QP);
  - the strand CLASS filters the arm-derived label: antiparallel ('+-'/'-+')
    drops cross-arm labels, congruent ('++'/'--') drops same-arm labels;
    contradictions demote to UNR with UNRESOLVED_TYPE=<label>_MISMATCH. This
    is a heuristic of resolve_single_tloc's own, not a rule
    resolve_translocation imposes (it only requires complementary strand
    classes across the two mates); see its docstring for the measurement;
  - CANDIDATE_SINGLE_TLOC is dispatched in the first pass of resolve() only:
    a cluster that shrinks to a single tloc after SR-only removal stays
    unresolved, identical to the flag-off outcome;
  - non-manta/dragen callers, missing EVIDENCE/PE and SR-only-only support
    stay unresolved exactly as before;
  - a breakpoint the cytoband table cannot answer for (CHR2 with no cytoband
    row; END2 past the last row) demotes to CTX_UNR instead of raising out of
    resolve, whatever exception pysam happens to raise;
  - the CLI wiring itself: --resolve-single-tlocs is accepted by argv and
    reaches ComplexSV, and a CTX record appears only with it.

Self-contained: generates its own fixtures in a temp dir, then exercises
ComplexSV directly (the same way svtk resolve does per cluster), the
_merge_records/sanity-check mechanics of the resolve CLI, and finally the CLI
end to end: svtk resolve's main() is run over the same fixture twice, with and
without --resolve-single-tlocs, so the wiring that carries the flag from argv
down to ComplexSV is covered, not just the flag's effect.

Requires: a python env with svtk importable (sys.path src/svtk) plus pysam,
pybedtools, and bgzip/tabix/bcftools on PATH.
"""
import contextlib
import io
import os
import shutil
import subprocess
import sys
import tempfile
from collections import deque
from importlib import util as importlib_util

import pysam
import pybedtools as pbt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'src', 'svtk'))
from svtk.cxsv import ComplexSV  # noqa: E402
from svtk.cxsv.complex_sv import get_arms  # noqa: E402

D = tempfile.mkdtemp(prefix='manta_tloc_test_')

# ---- fixtures ----
# cytobands (arm is field 4): chr1/chr2 with p/q boundaries at 50M/40M
with open(os.path.join(D, 'cytobands.bed'), 'w') as f:
    f.write('chr1\t0\t50000000\tp\n'
            'chr1\t50000000\t249250621\tq\n'
            'chr2\t0\t40000000\tp\n'
            'chr2\t40000000\t243199373\tq\n')
open(os.path.join(D, 'mei.bed'), 'w').close()
open(os.path.join(D, 'disc.bed'), 'w').close()

HEADER = [
    '##fileformat=VCFv4.2',
    '##source=svtk_test',
    '##contig=<ID=chr1,length=249250621>',
    '##contig=<ID=chr2,length=243199373>',
    # declared so that a record can name a CHR2 that has no cytoband row
    '##contig=<ID=chr9,length=138394717>',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant.">',
    '##INFO=<ID=CHR2,Number=1,Type=String,Description="Second position.">',
    '##INFO=<ID=END2,Number=1,Type=Integer,Description="Position of breakpoint on CHR2">',
    '##INFO=<ID=STRANDS,Number=1,Type=String,Description="Strand of the two breakpoints.">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the variant.">',
    '##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description="Algorithms which called the variant.">',
    '##INFO=<ID=EVIDENCE,Number=.,Type=String,Description="Classes of random forest support.">',
    '##INFO=<ID=MEMBERS,Number=.,Type=String,Description="IDs of constituent records.">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878']


def rec_line(vid, chrom, pos, alt, info):
    return f'{chrom}\t{pos}\t{vid}\tN\t{alt}\t.\t.\t{info}\tGT\t1/1'


# (vid, local pos, mate pos, alt, strands or None, algorithm, evidence or None)
DATA = [
    # --- resolving cases: strand class and arms agree ---
    # antiparallel +- , q<->q (same arm): paracentric -> CTX_PP/QQ
    ('tlocA', 100000000, 100000000, 'N[chr2:100000000[', '+-', 'manta', 'PE'),
    # antiparallel -+ , q<->q (same arm, other orientation) -> CTX_PP/QQ
    ('tlocB', 105000000, 105000000, ']chr2:105000000]N', '-+', 'manta', 'PE'),
    # congruent ++ , q<->p (different arm): pericentric -> CTX_PQ/QP
    # (this class was DROPPED by the first cut of the feature)
    ('tlocE', 106000000, 20000000, 'N]chr2:20000000]', '++', 'manta', 'PE'),
    # congruent -- , q<->p, dragen -> CTX_PQ/QP
    ('tlocM', 107000000, 21000000, ']chr2:21000000[N', '--', 'dragen', 'PE'),
    # antiparallel +- , q<->q, dragen -> CTX_PP/QQ (dragen gated like manta)
    ('tlocG', 110000000, 110000000, 'N[chr2:110000000[', '+-', 'dragen', 'PE'),
    # merged multi-caller site containing manta: gate is intersection-based,
    # documents that such sites DO resolve (the manta event is present)
    ('tlocP', 111000000, 111000000, 'N[chr2:111000000[', '+-', 'manta,wham', 'PE'),
    # --- demotions: arms contradict the strand class -> *_MISMATCH ---
    # antiparallel +- but cross-arm: paired path would not resolve either
    ('tlocC', 100000000, 20000000, 'N[chr2:20000000[', '+-', 'manta', 'PE'),
    # congruent ++ but same-arm
    ('tlocN', 120000000, 120000000, 'N]chr2:120000000]', '++', 'manta', 'PE'),
    # no STRANDS at all -> STRAND_MISMATCH_TLOC, no KeyError
    ('tlocO', 121000000, 121000000, 'N[chr2:121000000[', None, 'manta', 'PE'),
    # --- demotions mirroring the paired path, unchanged behavior ---
    # no EVIDENCE / SR-only EVIDENCE: CTX_UNR
    ('tlocI', 140000000, 140000000, 'N[chr2:140000000[', '+-', 'manta', None),
    ('tlocJ', 141000000, 141000000, 'N[chr2:141000000[', '+-', 'manta', 'SR'),
    # wham tloc: gated out by ALGORITHMS, stays UNR SINGLE_ENDER
    ('tlocD', 150000000, 150000000, 'N[chr2:150000000[', '+-', 'wham', 'PE'),
    # --- paired (two-mate) path must be untouched ---
    # genuine two-mate cluster at identical coords: DUPLICATE_COORDS guard
    ('tlocH1', 130000000, 130000000, 'N[chr2:130000000[', '+-', 'manta', 'PE'),
    ('tlocH2', 130000000, 130000000, ']chr2:130000000]N', '-+', 'manta', 'PE'),
    # manta PE tloc clustered with an SR-only wham ++ breakend: first pass
    # is STRAND_MISMATCH_TLOC; second pass removes the SR-only record and
    # must still auto-resolve the surviving single tloc
    ('tlocK', 160000000, 160000000, 'N[chr2:160000000[', '+-', 'manta', 'PE'),
    ('tlocL', 160000001, 160000000, 'N]chr2:160000000]', '++', 'wham', 'SR'),
]

lines = HEADER[:]
for vid, pos, mate, alt, strands, alg, ev in DATA:
    end2 = alt.split(':')[1].split(']')[0].split('[')[0]
    info = f'SVTYPE=BND;CHR2=chr2;END2={end2};SVLEN=-1;ALGORITHMS={alg};'
    if strands is not None:
        info += f'STRANDS={strands};'
    if ev is not None:
        info += f'EVIDENCE={ev};'
    info += f'MEMBERS={vid}'
    lines.append(rec_line(vid, 'chr1', pos, alt, info))
# simple insertion: unchanged path
lines.append(rec_line(
    'insF', 'chr1', 200000000, '<INS>',
    'SVTYPE=INS;CHR2=chr1;SVLEN=1000;STRANDS=+-;ALGORITHMS=manta;EVIDENCE=PE;MEMBERS=insF'))
# --- get_arms failures: the record is a well-formed manta tloc candidate, but
# the cytoband table cannot answer for one breakpoint. chr9 has no cytoband row
# at all; 300000000 is past the last chr2 row (chr2's q row ends 243199373).
# Both must demote to CTX_UNR, not raise out of resolve.
lines.append(rec_line(
    'tlocQ', 'chr1', 170000000, 'N[chr9:100000000[',
    'SVTYPE=BND;CHR2=chr9;END2=100000000;SVLEN=-1;ALGORITHMS=manta;STRANDS=+-;'
    'EVIDENCE=PE;MEMBERS=tlocQ'))
lines.append(rec_line(
    'tlocR', 'chr1', 171000000, 'N[chr2:300000000[',
    'SVTYPE=BND;CHR2=chr2;END2=300000000;SVLEN=-1;ALGORITHMS=manta;STRANDS=+-;'
    'EVIDENCE=PE;MEMBERS=tlocR'))

with open(os.path.join(D, 'raw.vcf'), 'w') as f:
    f.write('\n'.join(lines) + '\n')

for bed in ('cytobands', 'mei', 'disc'):
    subprocess.run(['bgzip', '-f', os.path.join(D, f'{bed}.bed')], check=True)
subprocess.run(['tabix', '-p', 'bed', os.path.join(D, 'cytobands.bed.gz')], check=True)
# svtk resolve opens --discfile with pysam.TabixFile, so it needs an index even
# though no case here triggers the single-INV rescan that would read it
subprocess.run(['tabix', '-p', 'bed', os.path.join(D, 'disc.bed.gz')], check=True)
subprocess.run(f'bgzip -c {os.path.join(D, "raw.vcf")} > {os.path.join(D, "raw.vcf.gz")}',
               shell=True, check=True)

# CLI-only copy: the same records minus tlocO. resolve.py's
# cluster_single_cleanup reads i.info['STRANDS'] unguarded, one step BEFORE
# ComplexSV is constructed, so a STRANDS-less record is a KeyError inside the
# CLI whatever the flag. tlocO exists to pin ComplexSV's own missing-STRANDS
# branch, so it stays out of the CLI run (and that resolve.py behaviour is a
# pre-existing property of main(), not of this feature).
with open(os.path.join(D, 'raw_cli.vcf'), 'w') as f:
    f.write('\n'.join(l for l in lines if '\ttlocO\t' not in l) + '\n')
subprocess.run(f'bgzip -c {os.path.join(D, "raw_cli.vcf")} > {os.path.join(D, "raw_cli.vcf.gz")}',
               shell=True, check=True)

# ---- load resolve CLI module (for CPX_INFO header lines & _merge_records) ----
spec = importlib_util.spec_from_file_location(
    'resolve_mod', os.path.join(HERE, 'src', 'svtk', 'svtk', 'cli', 'resolve.py'))
m = importlib_util.module_from_spec(spec)
spec.loader.exec_module(m)

cytobands = pysam.TabixFile(os.path.join(D, 'cytobands.bed.gz'))
mei_bed = pbt.BedTool(os.path.join(D, 'mei.bed.gz'))

src = pysam.VariantFile(os.path.join(D, 'raw.vcf.gz'))
for line in m.CPX_INFO:
    src.header.add_line(line)


def load_records():
    # fresh parse per case: ComplexSV stamps UNRESOLVED_TYPE etc. onto the
    # records it is given, so cases must not share mutable record objects
    vcf = pysam.VariantFile(os.path.join(D, 'raw.vcf.gz'))
    for line in m.CPX_INFO:
        vcf.header.add_line(line)
    return {r.id: r for r in vcf}


records = load_records()

# ---- cases ----
results = []


def run(vid, name, expected_svtype, expected_cpx, paired_with=None,
        expect_rec_unr=None, resolve_single_tlocs=True):
    recs = load_records()
    recs = [recs[vid]] if paired_with is None else [recs[vid], recs[paired_with]]
    cpx = ComplexSV(recs, cytobands, mei_bed, 1000,
                    resolve_single_tlocs=resolve_single_tlocs)
    ok = (cpx.svtype == expected_svtype) and (cpx.cpx_type == expected_cpx)
    extra = ''
    if cpx.svtype == 'CTX':
        v = cpx.vcf_record
        extra = (f' chrom={v.chrom}:{v.pos} CHR2={v.info.get("CHR2")} '
                 f'END2={v.info.get("END2")} ALT={v.alts}')
        ok = (ok and v.info.get('CHR2') == 'chr2' and v.alts == ('<CTX>',)
              and v.info.get('SVTYPE') == 'CTX'
              and v.info.get('CPX_TYPE') == expected_cpx
              and 'UNRESOLVED_TYPE' not in v.info.keys())
        if v.info['END2'] < v.pos:
            # production pysam 0.15.4 + htslib 1.9 write END=END2 (<POS) with
            # no clamp, identical to the paired path's stop=plus_end; pysam
            # >=0.22 silently drops stop<POS locally, so the assertion only
            # runs on the pinned stack
            if pysam.__version__.startswith('0.15'):
                ok = ok and v.stop == v.info['END2']
            else:
                print(f'NOTE {name}: local pysam {pysam.__version__} clamps '
                      f'stop<POS; production pysam 0.15.4 writes END=END2 '
                      f'(parity with paired path verified by code trace)')
    if expect_rec_unr is not None:
        stamped = cpx.records[0].info.get('UNRESOLVED_TYPE')
        ok = ok and stamped == expect_rec_unr
        extra += f' rec.UNRESOLVED_TYPE={stamped}'
    print(f'{"PASS" if ok else "FAIL"} {name}: svtype={cpx.svtype} '
          f'cpx_type={cpx.cpx_type}{extra}')
    results.append(ok)


# resolving: strand class consistent with arms
run('tlocA', 'manta +-  same-arm  ', 'CTX', 'CTX_PP/QQ')
run('tlocB', 'manta -+  same-arm  ', 'CTX', 'CTX_PP/QQ')
run('tlocE', 'manta ++  cross-arm ', 'CTX', 'CTX_PQ/QP')
run('tlocM', 'dragen -- cross-arm ', 'CTX', 'CTX_PQ/QP')
run('tlocG', 'dragen +- same-arm  ', 'CTX', 'CTX_PP/QQ')
run('tlocP', 'merged manta+wham   ', 'CTX', 'CTX_PP/QQ')
# opt-out (DEFAULT): every other svtk resolve consumer is unaffected -
# without the flag the identical record stays unresolved as in main
run('tlocA', 'flag OFF (default)  ', 'UNR', 'SINGLE_ENDER',
    resolve_single_tlocs=False)
# demotions: arms contradict strand class (record stamp must carry the
# suffixed type so the unresolved VCF is informative)
run('tlocC', 'manta +-  cross-arm ', 'UNR', 'CTX_PQ/QP_MISMATCH',
    expect_rec_unr='CTX_PQ/QP_MISMATCH')
run('tlocN', 'manta ++  same-arm  ', 'UNR', 'CTX_PP/QQ_MISMATCH',
    expect_rec_unr='CTX_PP/QQ_MISMATCH')
run('tlocO', 'manta strands absent', 'UNR', 'STRAND_MISMATCH_TLOC',
    expect_rec_unr='STRAND_MISMATCH_TLOC')
# demotions mirroring the paired path
run('tlocI', 'no EVIDENCE         ', 'UNR', 'CTX_UNR',
    expect_rec_unr='CTX_UNR')
run('tlocJ', 'SR-only EVIDENCE    ', 'UNR', 'CTX_UNR',
    expect_rec_unr='CTX_UNR')


# get_arms failure branch. Which exception pysam raises for an unanswerable
# breakpoint is not ours to assume: it depends on the shape of the bad
# coordinate AND on the pysam build (the image pins 0.15.4, this run may have a
# newer one). So measure it per case, then require that (i) the measured name
# is one of the three resolve_single_tloc actually catches, and (ii) ComplexSV
# absorbs it and demotes the record to CTX_UNR rather than raising.
# Measured here on pysam 0.24.1: a CHR2 with no cytoband row raises ValueError
# ('could not create iterator for region'), END2 past the last row raises
# StopIteration (the fetch is empty and get_arms takes next() of it). Narrowing
# the except tuple to any two of the three turns one of these two cases red;
# the third name, IndexError, is not reachable from either shape - it needs a
# cytoband row with fewer than four columns (get_arms does arm.split()[3][0]).
ARMS_CAUGHT = ('StopIteration', 'ValueError', 'IndexError')  # complex_sv.py


def arms_failure_case(vid, name, shape):
    recs = load_records()
    rec = recs[vid]
    try:
        arms = get_arms(rec, cytobands)
        raised = f'nothing (returned {arms})'
    except BaseException as exc:  # noqa: BLE001 - the type is the measurement
        arms, raised = None, type(exc).__name__
    escaped, cpx = None, None
    try:
        cpx = ComplexSV([rec], cytobands, mei_bed, 1000, resolve_single_tlocs=True)
    except BaseException as exc:  # noqa: BLE001 - must not escape resolve
        escaped = f'{type(exc).__name__}: {exc}'
    stamped = cpx.records[0].info.get('UNRESOLVED_TYPE') if cpx else None
    ok = (raised in ARMS_CAUGHT and escaped is None and cpx.svtype == 'UNR'
          and cpx.cpx_type == 'CTX_UNR' and stamped == 'CTX_UNR')
    print(f'{"PASS" if ok else "FAIL"} {name}: {shape} -> get_arms raised '
          f'{raised} (caught: {"/".join(ARMS_CAUGHT)}, pysam {pysam.__version__}); '
          f'svtype={cpx.svtype if cpx else None} cpx_type={cpx.cpx_type if cpx else None} '
          f'rec.UNRESOLVED_TYPE={stamped} escaped={escaped}')
    results.append(ok)


arms_failure_case('tlocQ', 'manta CHR2 no cytoband',
                  "CHR2=chr9, absent from the fixture cytoband table")
arms_failure_case('tlocR', 'manta END2 past last row',
                  'END2=300000000, past the last chr2 cytoband row (243199373)')

run('tlocD', 'wham tloc (gate)    ', 'UNR', 'SINGLE_ENDER')
run('insF', 'manta simple INS    ', 'INS', 'INS')
run('tlocH1', 'paired mates guard  ', 'UNR', 'CTX_PP/QQ_DUPLICATE_COORDS',
    paired_with='tlocH2')
# second-pass shrink: CANDIDATE_SINGLE_TLOC is dispatched in the FIRST pass
# only (svtk resolve hands the second pass INV single-enders only, never a
# BND). A cluster that only becomes a single tloc after the SR-only shrink is
# not re-dispatched: resolve() adds the SR-only record back and the cluster
# stays unresolved as STRAND_MISMATCH_TLOC, as it is with the flag off.
run('tlocK', 'second-pass shrink  ', 'UNR', 'STRAND_MISMATCH_TLOC',
    paired_with='tlocL', expect_rec_unr='STRAND_MISMATCH_TLOC')
# and the flag-on cluster must be field-for-field what flag-off produces for
# the same input, in either record order
for _order in (('tlocK', 'tlocL'), ('tlocL', 'tlocK')):
    def _state(flag, order=_order):
        c = ComplexSV([load_records()[k] for k in order], cytobands, mei_bed,
                      1000, resolve_single_tlocs=flag)
        v = c.vcf_record
        return (c.svtype, c.cpx_type, c.cluster_type, v.pos, v.alts,
                tuple(v.info.get('MEMBERS', ())),
                tuple(sorted(r.id for r in c.records)),
                tuple(r.info.get('UNRESOLVED_TYPE') for r in c.records))
    _on, _off = _state(True), _state(False)
    ok = _on == _off
    print(f'{"PASS" if ok else "FAIL"} shrink parity {_order} vs flag-off: '
          f'on={_on} off={_off}')
    results.append(ok)

# ---- CLI-flow simulation for the resolved CTX path in resolve.py main() ----
cpx = ComplexSV([load_records()['tlocE']], cytobands, mei_bed, 1000,
                resolve_single_tlocs=True)
cpx.vcf_record.id = 'CPX_CPX_chr1_000000'
cpx_record_ids = set(cpx.record_ids)
merged = list(m._merge_records(pysam.VariantFile(os.path.join(D, 'raw.vcf.gz')),
                               deque([cpx.vcf_record]), cpx_record_ids))
merged_ids = [r.id for r in merged]
ok = (merged_ids.count('CPX_CPX_chr1_000000') == 1 and 'tlocE' not in merged_ids
      and 'tlocE' in cpx.vcf_record.info['MEMBERS'])
print(f'{"PASS" if ok else "FAIL"} CLI merge: ctx emitted once, original '
      f'consumed, MEMBERS={cpx.vcf_record.info.get("MEMBERS")}')
results.append(ok)

used_vids = {r.id for r in merged}
for r in merged:
    used_vids.update(r.info['MEMBERS'] if 'MEMBERS' in r.info.keys() else ())
all_ids = {r.id for r in pysam.VariantFile(os.path.join(D, 'raw.vcf.gz'))}
ok = all_ids <= used_vids | cpx_record_ids
print(f'{"PASS" if ok else "FAIL"} CLI sanity check: consumed originals accounted for '
      f'(no POSTHOC_RESTORED duplicate of the resolved tloc)')
results.append(ok)

# ---- end-to-end CLI: the flag must reach ComplexSV through main() ----
# Everything above builds ComplexSV by hand, so the argv wiring in svtk
# resolve's main() is untested by them: delete the --resolve-single-tlocs
# argument (argparse then rejects the caller in mantatloccheck.sh) or delete
# resolve_single_tlocs=args.resolve_single_tlocs from the resolve_complex_sv
# call (production silently reverts to zero CTX), and every case above still
# prints ALL PASS. Run main() over this same fixture twice, with and without
# the flag, and compare the records it resolves.
ARMS_FAIL = {'tlocQ', 'tlocR'}


def run_resolve_cli(tag, extra_argv):
    """Run svtk resolve's main() in-process and capture what it resolved.

    The spy on resolve_complex_sv is what makes this assertable off the pinned
    image: with pysam >= 0.22 (0.24.1 here) main() opens its bcftools-sort
    writer from a header that has had CPX_INFO lines added but has not yielded
    a single record yet, and every write then dies with "Invalid BCF, the INFO
    tag id=16 is too large" - measured identically with the flag off, so it is
    a pysam-version artifact of running the CLI outside the image (which pins
    0.15.4), not a behaviour of this feature. The records main() pulled out of
    resolve_complex_sv are captured, so argv -> argparse -> resolve_complex_sv
    -> ComplexSV stays the real code path, and the output files are asserted
    too whenever the writes do succeed (the pinned stack).
    """
    captured = {'kwargs': None, 'records': []}
    real = m.resolve_complex_sv

    def spy(vcf, *args, **kwargs):
        captured['kwargs'] = kwargs
        for rec in real(vcf, *args, **kwargs):
            captured['records'].append(rec)
            yield rec

    resolved = os.path.join(D, f'{tag}.complex.vcf')
    unresolved = os.path.join(D, f'{tag}.unresolved.vcf')
    argv = [os.path.join(D, 'raw_cli.vcf.gz'), resolved,
            '--mei-bed', os.path.join(D, 'mei.bed.gz'),
            '--cytobands', os.path.join(D, 'cytobands.bed.gz'),
            '--discfile', os.path.join(D, 'disc.bed.gz'),
            '-u', unresolved] + extra_argv
    stopped = None
    m.resolve_complex_sv = spy
    try:
        with contextlib.redirect_stdout(io.StringIO()):  # main() is chatty
            m.main(argv)
    except BaseException as exc:  # noqa: BLE001 - SystemExit/write both results
        stopped = exc
    finally:
        m.resolve_complex_sv = real
    try:
        res = {r.id: r for r in pysam.VariantFile(resolved)}
    except Exception as exc:  # noqa: BLE001 - unpinned pysam wrote nothing usable
        res = None
    try:
        unr = {r.id: r.info.get('UNRESOLVED_TYPE') for r in pysam.VariantFile(unresolved)}
    except Exception:  # noqa: BLE001
        unr = None
    return captured, stopped, (res, unr)


def is_ctx(rec):
    return rec.alts == ('<CTX>',)


def outcome(records):
    """input record id -> what the CLI made of it. A resolved record is filed
    under the ids it consumed (MEMBERS, which for a resolved single record is
    that record's own id), because main() renames resolved records."""
    out = {}
    for r in records:
        if is_ctx(r):
            keys, state = r.info.get('MEMBERS', ()), 'CTX:' + str(r.info.get('CPX_TYPE'))
        elif r.info.get('UNRESOLVED'):
            keys, state = (r.id,), 'UNR:' + str(r.info.get('UNRESOLVED_TYPE'))
        else:
            keys = tuple(r.info.get('MEMBERS', ())) or (r.id,)
            state = 'RESOLVED:' + str(r.info.get('SVTYPE'))
        for k in keys:
            out[k] = state
    return out


run_on, stop_on, files_on = run_resolve_cli('flagon', ['--resolve-single-tlocs'])
run_off, stop_off, files_off = run_resolve_cli('flagoff', [])
seen_on, seen_off = outcome(run_on['records']), outcome(run_off['records'])

# 1. the flag survives argparse and reaches resolve_complex_sv both ways. A
# SystemExit here is argparse rejecting the caller's argv, i.e. exactly the
# regression this case exists for; anything else is the pysam write artifact
# described above, which stops main() after the records were already produced.
kw_on = run_on['kwargs'].get('resolve_single_tlocs') if run_on['kwargs'] else None
kw_off = run_off['kwargs'].get('resolve_single_tlocs') if run_off['kwargs'] else None
ok = (kw_on is True and not kw_off
      and not isinstance(stop_on, SystemExit) and not isinstance(stop_off, SystemExit)
      and run_on['records'] and run_off['records'])
print(f'{"PASS" if ok else "FAIL"} CLI flag wiring: resolve_single_tlocs='
      f'{kw_on} with the flag / {kw_off} without; '
      f'{len(run_on["records"])}/{len(run_off["records"])} records reached main(), '
      f'stopped={type(stop_on).__name__}/{type(stop_off).__name__}')
results.append(ok)

# 2. with the flag the CLI emits well-formed CTX records built out of the
# fixture's single manta tlocs. Which of them cluster together at CLI level is
# link_cpx's business (several DATA records sit at one locus on purpose, for the
# ComplexSV cases above), so the CTX set is read from this run rather than
# hard-coded; what is pinned is that it is non-empty, well-formed, drawn from
# the tlocs, and never built on an unanswerable breakpoint.
BND_IDS = {vid for vid, *_ in DATA} | ARMS_FAIL
ctx_on = [r for r in run_on['records'] if is_ctx(r)]
ctx_members = {m for r in ctx_on for m in r.info.get('MEMBERS', ())}
ok = (ctx_on
      and all(r.info['SVTYPE'] == 'CTX' and str(r.info.get('CPX_TYPE')).startswith('CTX_')
              and not r.info.get('UNRESOLVED') and r.info.get('CHR2') for r in ctx_on)
      and ctx_members <= BND_IDS and not ARMS_FAIL & ctx_members
      and seen_on.get('insF') == 'RESOLVED:INS')
print(f'{"PASS" if ok else "FAIL"} CLI e2e flag ON : {len(ctx_on)} CTX records '
      f'{sorted(ctx_members)}, insF={seen_on.get("insF")}')
results.append(ok)

# 3. the A/B: with the flag off there is no CTX at all, and every record the
# flag had consumed into a CTX comes back as an unresolved SINGLE_ENDER (not
# dropped). Same rule as case 4 of the paired-path cases: no END assertions.
ok = (not [r for r in run_off['records'] if is_ctx(r)]
      and ctx_members and not ctx_members.isdisjoint(set(seen_off))
      and all(str(seen_off.get(v, '')).startswith('UNR:') for v in ctx_members)
      and set(seen_off) == set(seen_on)
      and seen_off.get('insF') == 'RESOLVED:INS'
      and seen_off.get('tlocD') == seen_on.get('tlocD') == 'UNR:SINGLE_ENDER')
print(f'{"PASS" if ok else "FAIL"} CLI e2e flag OFF: 0 CTX, previously-consumed '
      f'{sorted(ctx_members)} now {[seen_off.get(v) for v in sorted(ctx_members)]}, '
      f'same record set={set(seen_off) == set(seen_on)}, wham/ins unchanged='
      f'{seen_off.get("tlocD")}/{seen_off.get("insF")}')
results.append(ok)

# 4. get_arms failures end to end: with the flag on the two unanswerable
# breakpoints do reach get_arms inside resolve and must come back as CTX_UNR
# unresolved records rather than raising out of the CLI (flag-off they never
# enter resolve_single_tloc, so they read as plain SINGLE_ENDER)
ok = (all(seen_on.get(v) == 'UNR:CTX_UNR' for v in ARMS_FAIL)
      and all(seen_off.get(v) == 'UNR:SINGLE_ENDER' for v in ARMS_FAIL)
      and not ARMS_FAIL & ctx_members)
print(f'{"PASS" if ok else "FAIL"} CLI e2e get_arms failure: flag-on '
      f'{[(v, seen_on.get(v)) for v in sorted(ARMS_FAIL)]}, flag-off '
      f'{[(v, seen_off.get(v)) for v in sorted(ARMS_FAIL)]}')
results.append(ok)

# 5. when the run really writes its two files (pinned pysam 0.15.4) check them
# too: CTX only in the flag-on resolved VCF
res_on, unr_on = files_on
res_off, unr_off = files_off
if res_on is None or unr_on is None or res_off is None or unr_off is None:
    print(f'NOTE CLI e2e output files unreadable on pysam {pysam.__version__}; '
          f'assertions 1-4 already cover argv -> ComplexSV. The pinned image '
          f'(pysam 0.15.4) writes both files and runs this check.')
else:
    ok = (sum(is_ctx(r) for r in res_on.values()) == len(ctx_on)
          and not any(is_ctx(r) for r in res_off.values())
          and all(unr_off.get(v) == 'SINGLE_ENDER' for v in ctx_members))
    print(f'{"PASS" if ok else "FAIL"} CLI e2e output files: CTX in resolved '
          f'flag-on={sum(is_ctx(r) for r in res_on.values())} '
          f'flag-off={sum(is_ctx(r) for r in res_off.values())}')
    results.append(ok)

shutil.rmtree(D, ignore_errors=True)
print()
if all(results):
    print('ALL PASS')
else:
    print('FAILURES PRESENT')
    sys.exit(1)
