package org.broadinstitute.svpipeline;

import java.io.ByteArrayInputStream;
import java.util.Arrays;

import org.broadinstitute.svpipeline.StitchFragmentedCNVs.PaddedInterval;
import org.broadinstitute.svpipeline.VCFParser.*;

public final class StitchFragmentedCNVsUnitTest {
    public static void main( final String[] args ) {
        testAsserts();
        testIsStitchable();
        testDoneStitching();
        testStitchTo();
        System.out.println("OK");
    }

    public static void testAsserts() {
        boolean caughtIt = false;
        try {
            assert (false);
        } catch ( final AssertionError ae ) {
            caughtIt = true;
        }
        if ( !caughtIt ) {
            throw new AssertionError("assertions aren't turned on (with -ea), so you're not testing anything.");
        }
    }

    public static void testIsStitchable() {
        final String vcfLine = "chr1\t1000\tID1\tN\t<DEL>\t60\tPASS\tEND=1999;SVTYPE=DEL;EVIDENCE=RD\tGT\t0/0\t0/1\n";
        final Record record = fromString(vcfLine);
        assert(StitchFragmentedCNVs.isStitchable(record));

        // not stitchable if there's a MULTIALLELIC filter component
        final ByteSequence originalFilter = record.getFilter().getValue();
        record.setFilter(StitchFragmentedCNVs.MULTIALLELIC);
        assert(!StitchFragmentedCNVs.isStitchable(record));
        final ByteSequence x = new ByteSequence("X");
        record.setFilter(Arrays.asList(x, StitchFragmentedCNVs.MULTIALLELIC));
        assert(!StitchFragmentedCNVs.isStitchable(record));
        final ByteSequence y = new ByteSequence("Y");
        record.setFilter(Arrays.asList(x, StitchFragmentedCNVs.MULTIALLELIC, y));
        assert(!StitchFragmentedCNVs.isStitchable(record));
        record.setFilter(originalFilter);
        assert(StitchFragmentedCNVs.isStitchable(record));

        // not stitchable if the the SVTYPE isn't DUP or DEL
        final InfoField info = record.getInfo();
        final ByteSequence originalSVTYPE = info.get(StitchFragmentedCNVs.SVTYPE);
        info.put(StitchFragmentedCNVs.SVTYPE, new ByteSequence("INS"));
        assert(!StitchFragmentedCNVs.isStitchable(record));
        info.put(StitchFragmentedCNVs.SVTYPE, StitchFragmentedCNVs.SVTYPE_DEL);
        assert(StitchFragmentedCNVs.isStitchable(record));
        info.put(StitchFragmentedCNVs.SVTYPE, StitchFragmentedCNVs.SVTYPE_DUP);
        assert(StitchFragmentedCNVs.isStitchable(record));
        info.put(StitchFragmentedCNVs.SVTYPE, originalSVTYPE);
        assert(StitchFragmentedCNVs.isStitchable(record));

        // not stitchable if END is missing
        final ByteSequence originalEnd = info.get(StitchFragmentedCNVs.END);
        info.remove(StitchFragmentedCNVs.END);
        assert(!StitchFragmentedCNVs.isStitchable(record));
        info.put(StitchFragmentedCNVs.END, originalEnd);
        assert(StitchFragmentedCNVs.isStitchable(record));

        // not stitchable if EVIDENCE includes PE or SR
        final ByteSequence originalEVIDENCE = info.get(StitchFragmentedCNVs.EVIDENCE);
        info.put(StitchFragmentedCNVs.EVIDENCE, StitchFragmentedCNVs.EVIDENCE_RD);
        assert(StitchFragmentedCNVs.isStitchable(record));
        info.put(StitchFragmentedCNVs.EVIDENCE, StitchFragmentedCNVs.EVIDENCE_BAF);
        assert(StitchFragmentedCNVs.isStitchable(record));
        info.put(StitchFragmentedCNVs.EVIDENCE, StitchFragmentedCNVs.EVIDENCE_SR);
        assert(!StitchFragmentedCNVs.isStitchable(record));
        info.put(StitchFragmentedCNVs.EVIDENCE, StitchFragmentedCNVs.EVIDENCE_PE);
        assert(!StitchFragmentedCNVs.isStitchable(record));
        final ByteSequence sep = new ByteSequence(",");
        info.put(StitchFragmentedCNVs.EVIDENCE,
                new ByteSequence(StitchFragmentedCNVs.EVIDENCE_BAF, sep, StitchFragmentedCNVs.EVIDENCE_SR));
        assert(!StitchFragmentedCNVs.isStitchable(record));
        info.put(StitchFragmentedCNVs.EVIDENCE, originalEVIDENCE);
        assert(StitchFragmentedCNVs.isStitchable(record));
    }

    public static void testDoneStitching() {
        final String vcfLine = "chr1\t1000\tID1\tN\t<DEL>\t60\tPASS\tEND=1999;SVTYPE=DEL;EVIDENCE=RD\tGT\t0/0\t0/1\n";
        final Record upstreamRecord = fromString(vcfLine);
        upstreamRecord.setPosition(1000);
        upstreamRecord.getInfo().put(StitchFragmentedCNVs.END, new ByteSequence(Integer.toString(1999)));
        final PaddedInterval upstreamInterval = new PaddedInterval(upstreamRecord);

        // we've got a record on chr1:1000-2000 with a 200-base pad.  move MAX_PAD bases further downstream.
        final int startNotTooFar = upstreamInterval.getPaddedEnd() + StitchFragmentedCNVs.MAX_PAD;
        final Record downstreamRecord = fromString(vcfLine);
        downstreamRecord.setPosition(startNotTooFar);
        downstreamRecord.getInfo().put(StitchFragmentedCNVs.END, new ByteSequence(Integer.toString(startNotTooFar + 1000)));
        final PaddedInterval downstreamInterval = new PaddedInterval(downstreamRecord);
        assert(!upstreamInterval.doneStitching(downstreamInterval));

        // move one more base, and doneStitching should return true
        downstreamRecord.setPosition(startNotTooFar + 1);
        assert(upstreamInterval.doneStitching(new PaddedInterval(downstreamRecord)));
    }

    public static void testStitchTo() {
        final String vcfLine1 = "chr1\t1000\tID1\tN\t<DEL>\t60\tPASS\tEND=1999;SVTYPE=DEL;EVIDENCE=RD\tGT\t0/0\t0/1\n";
        final Record upstreamRecord = fromString(vcfLine1);
        assert(StitchFragmentedCNVs.isStitchable(upstreamRecord));
        final PaddedInterval upstreamInterval = new PaddedInterval(upstreamRecord);
        final String vcfLine2 = "chr1\t2399\tID1\tN\t<DEL>\t60\tPASS\tEND=3399;SVTYPE=DEL;EVIDENCE=RD\tGT\t0/0\t0/1\n";
        final Record downstreamRecord = fromString(vcfLine2);
        assert(StitchFragmentedCNVs.isStitchable(downstreamRecord));
        final PaddedInterval stitched = upstreamInterval.stitchTo(new PaddedInterval(downstreamRecord));
        assert(stitched != null);
        assert(stitched.getRecord().getPosition() == 1000);
        assert(stitched.getVCFEnd() == 3399);

        // fails because no overlap (padded intervals are adjacent)
        downstreamRecord.setPosition(2400);
        assert(upstreamInterval.stitchTo(new PaddedInterval(downstreamRecord)) == null);

        // back to starting conditions
        downstreamRecord.setPosition(2399);
        assert(upstreamInterval.stitchTo(new PaddedInterval(downstreamRecord)) != null);

        // fails because event types don't match
        downstreamRecord.getInfo().put(StitchFragmentedCNVs.SVTYPE, StitchFragmentedCNVs.SVTYPE_DUP);
        assert(upstreamInterval.stitchTo(new PaddedInterval(downstreamRecord)) == null);

        // back to starting conditions
        downstreamRecord.getInfo().put(StitchFragmentedCNVs.SVTYPE, StitchFragmentedCNVs.SVTYPE_DEL);
        assert(upstreamInterval.stitchTo(new PaddedInterval(downstreamRecord)) != null);

        // overlaps upstream interval too much
        downstreamRecord.setPosition(1799);
        assert(upstreamInterval.stitchTo(new PaddedInterval(downstreamRecord)) == null);

        // back to starting conditions
        downstreamRecord.setPosition(2399);
        assert(upstreamInterval.stitchTo(new PaddedInterval(downstreamRecord)) != null);

        // overlaps downstream interval too much
        downstreamRecord.setPosition(1899);
        downstreamRecord.getInfo().put(StitchFragmentedCNVs.END, new ByteSequence(Integer.toString(2399)));
        assert(upstreamInterval.stitchTo(new PaddedInterval(downstreamRecord)) == null);

        // back to starting conditions
        downstreamRecord.setPosition(2399);
        downstreamRecord.getInfo().put(StitchFragmentedCNVs.END, new ByteSequence(Integer.toString(3399)));
        assert(upstreamInterval.stitchTo(new PaddedInterval(downstreamRecord)) != null);

        // genotypes don't match
        downstreamRecord.getGenotypes().get(0).set(0, new ByteSequence("0/1"));
        assert(upstreamInterval.stitchTo(new PaddedInterval(downstreamRecord)) == null);
    }

    private static Record fromString( final String vcfLine ) {
        return new VCFParser(new ByteArrayInputStream(vcfLine.getBytes())).nextRecord();
    }
}
