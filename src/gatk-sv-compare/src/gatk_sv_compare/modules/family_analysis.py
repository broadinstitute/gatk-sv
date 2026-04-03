"""Family-based inheritance and de novo-rate diagnostics."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from ..dimensions import ordered_plot_af_buckets, ordered_plot_size_buckets, ordered_svtypes
from ..plot_utils import SUMMARY_COLORS, SVTYPE_COLORS, double_column_figsize, single_column_figsize, plot_beeswarm_horizontal, plot_heatmap_annotated, save_figure
from .base import AnalysisModule, write_tsv_gz


@dataclass(frozen=True)
class Trio:
    family_id: str
    proband: str
    father: Optional[str]
    mother: Optional[str]

    @property
    def family_type(self) -> str:
        return "trios" if self.father and self.mother else "duos"


@dataclass(frozen=True)
class TransmissionRecord:
    label: str
    family_id: str
    family_type: str
    variant_id: str
    svtype: str
    size_bucket: str
    af_bucket: str
    contig: str
    proband_alt: int
    proband_inherited: float
    proband_denovo: float
    father_alt: int
    father_transmitted: float
    father_untransmitted: float
    mother_alt: int
    mother_transmitted: float
    mother_untransmitted: float
    proband_gq: Optional[float]
    father_gq: Optional[float]
    mother_gq: Optional[float]


def parse_ped_file(ped_path: Path, vcf_samples: set[str]) -> List[Trio]:
    trios: List[Trio] = []
    for line in ped_path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) < 4:
            continue
        family_id, proband, father, mother = fields[:4]
        if proband not in vcf_samples:
            continue
        father_id = None if father in {"0", ".", "NA"} or father not in vcf_samples else father
        mother_id = None if mother in {"0", ".", "NA"} or mother not in vcf_samples else mother
        if father_id is None and mother_id is None:
            continue
        trios.append(Trio(family_id=family_id, proband=proband, father=father_id, mother=mother_id))
    return trios


def _record_identifier(record: pysam.VariantRecord) -> str:
    return record.id or f"{record.contig}:{record.pos}"


def _iter_records_for_contig(vcf: pysam.VariantFile, contig: str):
    try:
        yield from vcf.fetch(contig)
    except ValueError:
        for record in vcf:
            if record.contig == contig:
                yield record


def _is_autosome(contig: str) -> bool:
    normalized = contig[3:] if contig.startswith("chr") else contig
    return normalized.isdigit() and 1 <= int(normalized) <= 22


def _gt_alt_count(sample: pysam.libcbcf.VariantRecordSample) -> Optional[int]:
    gt = sample.get("GT")
    if not gt or any(allele is None for allele in gt):
        return None
    return int(sum(1 for allele in gt if allele and allele > 0))


def _gq(sample: pysam.libcbcf.VariantRecordSample) -> Optional[float]:
    value = sample.get("GQ")
    if value in (None, "."):
        return None
    return float(value)


def _transmission_from_counts(proband: int, father: int, mother: int) -> Dict[str, float]:
    parental_total = father + mother
    proband_denovo = float(max(0, proband - parental_total))
    proband_inherited = float(max(0, proband - proband_denovo))
    if parental_total > 0 and proband_inherited > 0:
        father_transmitted = proband_inherited * (father / parental_total)
        mother_transmitted = proband_inherited * (mother / parental_total)
    else:
        father_transmitted = 0.0
        mother_transmitted = 0.0
    return {
        "proband_inherited": proband_inherited,
        "proband_denovo": proband_denovo,
        "father_transmitted": father_transmitted,
        "father_untransmitted": float(max(0, father - father_transmitted)),
        "mother_transmitted": mother_transmitted,
        "mother_untransmitted": float(max(0, mother - mother_transmitted)),
    }


def build_transmission_table(
    vcf_path: Path,
    sites: pd.DataFrame,
    label: str,
    trios: Sequence[Trio],
    pass_only: bool = False,
) -> pd.DataFrame:
    if not trios:
        return pd.DataFrame(columns=[field for field in TransmissionRecord.__dataclass_fields__.keys()])
    filtered_sites = sites.copy()
    filtered_sites = filtered_sites.loc[filtered_sites["svtype"].isin(["DEL", "DUP", "INS", "INS:MEI", "INV"])]
    filtered_sites = filtered_sites.loc[filtered_sites["contig"].map(_is_autosome)]
    if pass_only:
        filtered_sites = filtered_sites.loc[filtered_sites["in_filtered_pass_view"]]
    if filtered_sites.empty:
        return pd.DataFrame(columns=[field for field in TransmissionRecord.__dataclass_fields__.keys()])

    site_meta = filtered_sites.set_index("variant_id")[["svtype", "size_bucket", "af_bucket", "contig"]]
    target_ids_by_contig: Dict[str, set[str]] = {}
    for contig, frame in filtered_sites.groupby("contig"):
        target_ids_by_contig[str(contig)] = set(frame["variant_id"].astype(str))

    rows: List[dict] = []
    with pysam.VariantFile(str(vcf_path)) as vcf:
        sample_names = list(vcf.header.samples)
        sample_to_index = {sample: idx for idx, sample in enumerate(sample_names)}
        valid_families = [
            trio for trio in trios if trio.proband in sample_to_index and (trio.father is None or trio.father in sample_to_index) and (trio.mother is None or trio.mother in sample_to_index)
        ]
        for contig, target_ids in target_ids_by_contig.items():
            for record in _iter_records_for_contig(vcf, contig):
                variant_id = _record_identifier(record)
                if variant_id not in target_ids:
                    continue
                meta = site_meta.loc[variant_id]
                samples = list(record.samples.values())
                for trio in valid_families:
                    proband_sample = samples[sample_to_index[trio.proband]]
                    father_sample = samples[sample_to_index[trio.father]] if trio.father else None
                    mother_sample = samples[sample_to_index[trio.mother]] if trio.mother else None
                    proband_alt = _gt_alt_count(proband_sample)
                    father_alt = 0 if father_sample is None else _gt_alt_count(father_sample)
                    mother_alt = 0 if mother_sample is None else _gt_alt_count(mother_sample)
                    if proband_alt is None or (father_sample is not None and father_alt is None) or (mother_sample is not None and mother_alt is None):
                        continue
                    transmission = _transmission_from_counts(proband_alt, int(father_alt or 0), int(mother_alt or 0))
                    rows.append(
                        TransmissionRecord(
                            label=label,
                            family_id=trio.family_id,
                            family_type=trio.family_type,
                            variant_id=variant_id,
                            svtype=str(meta["svtype"]),
                            size_bucket=str(meta["size_bucket"]),
                            af_bucket=str(meta["af_bucket"]),
                            contig=str(meta["contig"]),
                            proband_alt=int(proband_alt),
                            proband_inherited=float(transmission["proband_inherited"]),
                            proband_denovo=float(transmission["proband_denovo"]),
                            father_alt=int(father_alt or 0),
                            father_transmitted=float(transmission["father_transmitted"]),
                            father_untransmitted=float(transmission["father_untransmitted"]),
                            mother_alt=int(mother_alt or 0),
                            mother_transmitted=float(transmission["mother_transmitted"]),
                            mother_untransmitted=float(transmission["mother_untransmitted"]),
                            proband_gq=_gq(proband_sample),
                            father_gq=None if father_sample is None else _gq(father_sample),
                            mother_gq=None if mother_sample is None else _gq(mother_sample),
                        ).__dict__
                    )
    return pd.DataFrame(rows)


def summarize_inheritance_stats(records: pd.DataFrame) -> pd.DataFrame:
    if records.empty:
        return pd.DataFrame(columns=["label", "family_type", "family_id", "site_proband_total", "site_proband_inherited", "site_proband_denovo", "site_de_novo_rate", "allele_proband_total", "allele_proband_inherited", "allele_proband_denovo", "allele_de_novo_rate", "site_paternal_fraction", "site_maternal_fraction", "allele_paternal_fraction", "allele_maternal_fraction"])

    rows = []
    for (label, family_type, family_id), group in records.groupby(["label", "family_type", "family_id"], dropna=False):
        site_proband_total = int((group["proband_alt"] > 0).sum())
        site_proband_inherited = int(((group["proband_alt"] > 0) & (group["proband_inherited"] > 0)).sum())
        site_proband_denovo = int(((group["proband_alt"] > 0) & (group["proband_denovo"] > 0)).sum())
        allele_proband_total = float(group["proband_alt"].sum())
        allele_proband_inherited = float(group["proband_inherited"].sum())
        allele_proband_denovo = float(group["proband_denovo"].sum())
        site_pat = float(group["father_transmitted"].gt(0).sum())
        site_mat = float(group["mother_transmitted"].gt(0).sum())
        rows.append(
            {
                "label": label,
                "family_type": family_type,
                "family_id": family_id,
                "site_proband_total": site_proband_total,
                "site_proband_inherited": site_proband_inherited,
                "site_proband_denovo": site_proband_denovo,
                "site_inheritance_rate": site_proband_inherited / site_proband_total if site_proband_total else np.nan,
                "site_de_novo_rate": site_proband_denovo / site_proband_total if site_proband_total else np.nan,
                "site_father_transmission_rate": site_pat / max(float(group["father_alt"].gt(0).sum()), 1.0),
                "site_mother_transmission_rate": site_mat / max(float(group["mother_alt"].gt(0).sum()), 1.0),
                "site_paternal_fraction": site_pat / max(site_pat + site_mat, 1.0),
                "site_maternal_fraction": site_mat / max(site_pat + site_mat, 1.0),
                "allele_proband_total": allele_proband_total,
                "allele_proband_inherited": allele_proband_inherited,
                "allele_proband_denovo": allele_proband_denovo,
                "allele_inheritance_rate": allele_proband_inherited / allele_proband_total if allele_proband_total else np.nan,
                "allele_de_novo_rate": allele_proband_denovo / allele_proband_total if allele_proband_total else np.nan,
                "allele_father_transmission_rate": float(group["father_transmitted"].sum()) / max(float(group["father_alt"].sum()), 1.0),
                "allele_mother_transmission_rate": float(group["mother_transmitted"].sum()) / max(float(group["mother_alt"].sum()), 1.0),
                "allele_paternal_fraction": float(group["father_transmitted"].sum()) / max(float(group["proband_inherited"].sum()), 1.0),
                "allele_maternal_fraction": float(group["mother_transmitted"].sum()) / max(float(group["proband_inherited"].sum()), 1.0),
            }
        )
    return pd.DataFrame(rows)


def summarize_denovo_by_dimension(records: pd.DataFrame, dimension: str) -> pd.DataFrame:
    if records.empty:
        return pd.DataFrame(columns=["label", "family_type", dimension, "svtype", "site_denovo_rate", "allele_denovo_rate", "n_records"])
    rows = []
    for keys, group in records.groupby(["label", "family_type", dimension, "svtype"], dropna=False):
        label, family_type, dim_value, svtype = keys
        site_total = float((group["proband_alt"] > 0).sum())
        site_denovo = float(((group["proband_alt"] > 0) & (group["proband_denovo"] > 0)).sum())
        allele_total = float(group["proband_alt"].sum())
        allele_denovo = float(group["proband_denovo"].sum())
        rows.append(
            {
                "label": label,
                "family_type": family_type,
                dimension: dim_value,
                "svtype": svtype,
                "site_denovo_rate": site_denovo / site_total if site_total else np.nan,
                "allele_denovo_rate": allele_denovo / allele_total if allele_total else np.nan,
                "n_records": len(group),
            }
        )
    return pd.DataFrame(rows)


def summarize_denovo_by_gq(records: pd.DataFrame, thresholds: Sequence[int]) -> pd.DataFrame:
    rows = []
    for threshold in thresholds:
        filtered = records.loc[records["proband_gq"].fillna(-1) >= float(threshold)]
        if filtered.empty:
            continue
        for (label, family_type, svtype), group in filtered.groupby(["label", "family_type", "svtype"], dropna=False):
            site_total = float((group["proband_alt"] > 0).sum())
            allele_total = float(group["proband_alt"].sum())
            rows.append(
                {
                    "label": label,
                    "family_type": family_type,
                    "svtype": svtype,
                    "min_gq": int(threshold),
                    "site_denovo_rate": float((((group["proband_alt"] > 0) & (group["proband_denovo"] > 0)).sum()) / site_total) if site_total else np.nan,
                    "allele_denovo_rate": float(group["proband_denovo"].sum() / allele_total) if allele_total else np.nan,
                    "n_records": len(group),
                }
            )
    return pd.DataFrame(rows)


def summarize_size_x_freq(records: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    outputs: Dict[str, pd.DataFrame] = {}
    for svtype, group in [("all", records)] + [(svtype, records.loc[records["svtype"] == svtype]) for svtype in ordered_svtypes(records["svtype"].unique())]:
        if group.empty:
            outputs[svtype] = pd.DataFrame(columns=["size_bucket", "af_bucket", "site_denovo_rate", "allele_denovo_rate", "n_records"])
            continue
        table = group.groupby(["size_bucket", "af_bucket"], dropna=False).agg(
            n_records=("variant_id", "count"),
            site_denovo_rate=("proband_denovo", lambda values: float((values > 0).mean())),
            allele_denovo_rate=("proband_denovo", lambda values: float(values.sum() / max(group.loc[values.index, "proband_alt"].sum(), 1.0))),
        ).reset_index()
        outputs[svtype] = table
    return outputs


def _plot_inheritance_beeswarm(stats: pd.DataFrame, family_type: str, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=double_column_figsize(3.0))
    subset = stats.loc[stats["family_type"] == family_type]
    if subset.empty:
        ax.text(0.5, 0.5, "No family data", ha="center", va="center")
        ax.set_axis_off()
    else:
        values_list = [
            subset["site_inheritance_rate"].dropna().tolist(),
            subset["site_de_novo_rate"].dropna().tolist(),
            subset["allele_inheritance_rate"].dropna().tolist(),
            subset["allele_de_novo_rate"].dropna().tolist(),
        ]
        colors = [SUMMARY_COLORS["tertiary"], SUMMARY_COLORS["quaternary"], SUMMARY_COLORS["primary"], SUMMARY_COLORS["secondary"]]
        labels = ["site inherited", "site de novo", "allele inherited", "allele de novo"]
        plot_beeswarm_horizontal(ax, values_list, colors, labels)
        ax.set_xlim(0.0, 1.0)
        ax.set_xlabel("Rate")
        ax.set_title(f"Inheritance summary: {family_type}")
    save_figure(fig, output_path)


def _plot_dnr_curve(table: pd.DataFrame, x_field: str, y_field: str, family_type: str, output_path: Path, log_x: bool) -> None:
    fig, ax = plt.subplots(figsize=double_column_figsize(3.0))
    subset = table.loc[table["family_type"] == family_type].copy()
    if subset.empty:
        ax.text(0.5, 0.5, "No family data", ha="center", va="center")
        ax.set_axis_off()
        save_figure(fig, output_path)
        return
    if x_field == "size_bucket":
        x_values = ordered_plot_size_buckets(subset[x_field].tolist())
    elif x_field == "af_bucket":
        x_values = ordered_plot_af_buckets(subset[x_field].tolist())
    else:
        x_values = list(dict.fromkeys(subset.sort_values(x_field)[x_field].tolist()))
    if all(isinstance(value, (int, float, np.number)) for value in x_values):
        plot_positions = [float(value) for value in x_values]
        if log_x and all(float(value) > 0 for value in x_values):
            ax.set_xscale("log")
    else:
        plot_positions = list(np.arange(len(x_values), dtype=float))
        ax.set_xticks(plot_positions)
        ax.set_xticklabels([str(value) for value in x_values], rotation=30)
        ax.set_xscale("linear")

    for svtype in ordered_svtypes(subset["svtype"].unique()):
        group = subset.loc[subset["svtype"] == svtype]
        ordered = group.groupby(x_field, dropna=False)[y_field].median()
        series = [float(ordered.get(value, np.nan)) for value in x_values]
        ax.plot(plot_positions, series, label=str(svtype), color=SVTYPE_COLORS.get(str(svtype), SUMMARY_COLORS["neutral"]))
    ax.legend()
    ax.set_ylabel(y_field.replace("_", " "))
    ax.set_title(f"{y_field.replace('_', ' ')} vs {x_field.replace('_', ' ')}: {family_type}")
    save_figure(fig, output_path)


def _plot_dnr_heatmap(table: pd.DataFrame, output_path: Path, title: str, value_field: str) -> None:
    fig, ax = plt.subplots(figsize=single_column_figsize(3.2))
    if table.empty:
        ax.text(0.5, 0.5, "No family data", ha="center", va="center")
        ax.set_axis_off()
        save_figure(fig, output_path)
        return
    matrix = table.pivot(index="size_bucket", columns="af_bucket", values=value_field).fillna(0.0)
    matrix = matrix.reindex(index=ordered_plot_size_buckets(matrix.index), columns=ordered_plot_af_buckets(matrix.columns), fill_value=0.0)
    image = plot_heatmap_annotated(
        ax,
        matrix.values,
        list(matrix.index),
        list(matrix.columns),
        fmt="{value:.2f}",
        value_range=(0.0, 1.0),
    )
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("De novo rate")
    ax.set_title(title)
    save_figure(fig, output_path)


class FamilyAnalysisModule(AnalysisModule):
    @property
    def name(self) -> str:
        return "family_analysis"

    @property
    def requires_ped_file(self) -> bool:
        return True

    @property
    def requires_genotype_pass(self) -> bool:
        return True

    def run(self, data: AggregatedData, config: AnalysisConfig) -> None:
        if config.ped_file is None:
            return
        trios = parse_ped_file(config.ped_file, set(data.sample_names_a) | set(data.sample_names_b))
        if not trios:
            return

        output_dir = self.output_dir(config)
        tables_dir = output_dir / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)

        records_a = build_transmission_table(config.vcf_a_path, data.sites_a, data.label_a, trios, pass_only=config.pass_only)
        records_b = build_transmission_table(config.vcf_b_path, data.sites_b, data.label_b, trios, pass_only=config.pass_only)
        records = pd.concat([records_a, records_b], ignore_index=True) if not records_a.empty or not records_b.empty else pd.DataFrame()
        if records.empty:
            return

        inheritance = summarize_inheritance_stats(records)
        write_tsv_gz(inheritance, tables_dir / "inheritance_stats.trios.tsv")

        by_class = summarize_denovo_by_dimension(records, "svtype")
        by_size = summarize_denovo_by_dimension(records, "size_bucket")
        by_freq = summarize_denovo_by_dimension(records, "af_bucket")
        max_gq = int(records["proband_gq"].fillna(0).max()) if not records.empty else 0
        thresholds = np.unique(np.linspace(0, max(max_gq, 1), num=min(max(max_gq, 1), 10), dtype=int))
        by_gq = summarize_denovo_by_gq(records, thresholds)
        write_tsv_gz(by_class, tables_dir / "denovo_rate_by_class.tsv")
        write_tsv_gz(by_size, tables_dir / "denovo_rate_by_size.tsv")
        write_tsv_gz(by_freq, tables_dir / "denovo_rate_by_freq.tsv")
        write_tsv_gz(by_gq, tables_dir / "denovo_rate_by_gq.tsv")

        for svtype, table in summarize_size_x_freq(records).items():
            write_tsv_gz(table, tables_dir / f"denovo_rate_size_x_freq.{svtype}.tsv")

        for family_type in sorted(records["family_type"].unique()):
            _plot_inheritance_beeswarm(inheritance, family_type, output_dir / f"inheritance.{family_type}.all_sv.png")
            _plot_dnr_curve(by_size, "size_bucket", "site_denovo_rate", family_type, output_dir / f"dnr_vs_size.{family_type}.variants.png", log_x=False)
            _plot_dnr_curve(by_size, "size_bucket", "allele_denovo_rate", family_type, output_dir / f"dnr_vs_size.{family_type}.alleles.png", log_x=False)
            _plot_dnr_curve(by_freq, "af_bucket", "site_denovo_rate", family_type, output_dir / f"dnr_vs_freq.{family_type}.variants.png", log_x=False)
            _plot_dnr_curve(by_freq, "af_bucket", "allele_denovo_rate", family_type, output_dir / f"dnr_vs_freq.{family_type}.alleles.png", log_x=False)
            _plot_dnr_curve(by_gq, "min_gq", "site_denovo_rate", family_type, output_dir / f"dnr_vs_gq.{family_type}.variants.png", log_x=False)
            _plot_dnr_curve(by_gq, "min_gq", "allele_denovo_rate", family_type, output_dir / f"dnr_vs_gq.{family_type}.alleles.png", log_x=False)
            overall_heat = summarize_size_x_freq(records).get("all", pd.DataFrame())
            _plot_dnr_heatmap(overall_heat, output_dir / f"dnr_heatmap.{family_type}.variants.all_sv.png", f"De novo heatmap: {family_type}", "site_denovo_rate")
