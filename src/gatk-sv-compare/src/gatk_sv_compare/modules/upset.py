"""UpSet plots for ALGORITHMS and EVIDENCE combinations."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from ..dimensions import EVIDENCE_COMPONENT_ORDER, normalize_algorithms, normalize_evidence_bucket, ordered_algorithms
from ..plot_utils import double_column_figsize, save_figure
from .base import AnalysisModule, write_tsv_gz

try:
    from upsetplot import UpSet, from_memberships
except ImportError:  # pragma: no cover - exercised only when dependency is absent
    UpSet = None
    from_memberships = None


if UpSet is not None:
    class CompatibleUpSet(UpSet):
        def plot_matrix(self, ax):
            ax = self._reorient(ax)
            data = self.intersections
            n_cats = data.index.nlevels

            inclusion = data.index.to_frame().values

            styles = [
                [
                    self.subset_styles[i]
                    if inclusion[i, j]
                    else {"facecolor": self._other_dots_color, "linewidth": 0}
                    for j in range(n_cats)
                ]
                for i in range(len(data))
            ]
            styles = sum(styles, [])
            style_columns = {
                "facecolor": "facecolors",
                "edgecolor": "edgecolors",
                "linewidth": "linewidths",
                "linestyle": "linestyles",
                "hatch": "hatch",
            }
            styles = pd.DataFrame(styles).reindex(columns=style_columns.keys()).astype(
                {
                    "facecolor": "O",
                    "edgecolor": "O",
                    "linewidth": float,
                    "linestyle": "O",
                    "hatch": "O",
                }
            )
            styles["linewidth"] = styles["linewidth"].fillna(1.0)
            styles["facecolor"] = styles["facecolor"].fillna(self._facecolor)
            styles["edgecolor"] = styles["edgecolor"].fillna(styles["facecolor"])
            styles["linestyle"] = styles["linestyle"].fillna("solid")
            del styles["hatch"]

            x = np.repeat(np.arange(len(data)), n_cats)
            y = np.tile(np.arange(n_cats), len(data))

            if self._element_size is not None:
                s = (self._element_size * 0.35) ** 2
            else:
                s = 200
            ax.scatter(
                *self._swapaxes(x, y),
                s=s,
                zorder=10,
                **styles.rename(columns=style_columns),
            )

            if self._with_lines:
                idx = np.flatnonzero(inclusion)
                line_data = pd.Series(y[idx], index=x[idx]).groupby(level=0).aggregate(["min", "max"])
                colors = pd.Series(
                    [
                        style.get("edgecolor", style.get("facecolor", self._facecolor))
                        for style in self.subset_styles
                    ],
                    name="color",
                )
                line_data = line_data.join(colors)
                ax.vlines(
                    line_data.index.values,
                    line_data["min"],
                    line_data["max"],
                    lw=2,
                    colors=line_data["color"],
                    zorder=5,
                )

            tick_axis = ax.yaxis
            tick_axis.set_ticks(np.arange(n_cats))
            tick_axis.set_ticklabels(data.index.names, rotation=0 if self._horizontal else -90)
            ax.xaxis.set_visible(False)
            ax.tick_params(axis="both", which="both", length=0)
            if not self._horizontal:
                ax.yaxis.set_ticks_position("top")
            ax.set_frame_on(False)
            ax.set_xlim(-0.5, x[-1] + 0.5, auto=False)
            ax.grid(False)
else:  # pragma: no cover - exercised only when dependency is absent
    CompatibleUpSet = None


@dataclass(frozen=True)
class MembershipField:
    name: str
    source_column: str
    title: str


_MEMBERSHIP_FIELDS = (
    MembershipField(name="algorithms", source_column="algorithms", title="Algorithms"),
    MembershipField(name="evidence", source_column="evidence_bucket", title="Evidence"),
)


def _require_upsetplot() -> None:
    if UpSet is None or from_memberships is None:
        raise RuntimeError("upsetplot is required to run the upset module. Install it with `pip install upsetplot`.")


def _filtered_sites(sites: pd.DataFrame, pass_only: bool) -> pd.DataFrame:
    if not pass_only:
        return sites.copy()
    return sites.loc[sites["in_filtered_pass_view"]].copy()


def _parse_membership(value: object, field: MembershipField) -> tuple[str, ...]:
    if field.name == "algorithms":
        return tuple(normalize_algorithms(value))
    evidence_bucket = normalize_evidence_bucket(value)
    if evidence_bucket == "unknown":
        return ("unknown",)
    return tuple(item for item in evidence_bucket.split(",") if item)


def _membership_category_order(memberships: Iterable[Sequence[str]], field: MembershipField) -> list[str]:
    observed = {item for membership in memberships for item in membership}
    if field.name == "algorithms":
        return [item for item in ordered_algorithms(observed) if item in observed]
    ordered = [item for item in EVIDENCE_COMPONENT_ORDER if item in observed]
    if "unknown" in observed:
        ordered.append("unknown")
    extras = [item for item in sorted(observed) if item not in set(ordered)]
    return ordered + extras


def build_membership_count_table(sites: pd.DataFrame, field_name: str, pass_only: bool = False) -> pd.DataFrame:
    field = next(config for config in _MEMBERSHIP_FIELDS if config.name == field_name)
    filtered = _filtered_sites(sites, pass_only)
    if filtered.empty:
        return pd.DataFrame(columns=["combination", "n_variants", "n_categories"])
    counts = Counter(_parse_membership(value, field) for value in filtered[field.source_column])
    rows = [
        {
            "combination": ",".join(membership),
            "n_variants": count,
            "n_categories": len(membership),
        }
        for membership, count in counts.items()
    ]
    result = pd.DataFrame(rows)
    return result.sort_values(["n_variants", "n_categories", "combination"], ascending=[False, False, True]).reset_index(drop=True)


def _build_upset_series(sites: pd.DataFrame, field: MembershipField, pass_only: bool = False) -> pd.Series | None:
    _require_upsetplot()
    filtered = _filtered_sites(sites, pass_only)
    if filtered.empty:
        return None
    counts = Counter(_parse_membership(value, field) for value in filtered[field.source_column])
    memberships = [list(membership) for membership in counts]
    values = [counts[tuple(membership)] for membership in counts]
    series = from_memberships(memberships, data=values)
    desired_order = [item for item in _membership_category_order(counts.keys(), field) if item in series.index.names]
    if desired_order and list(series.index.names) != desired_order:
        series = series.reorder_levels(desired_order)
    return series.sort_values(ascending=False)


def _plot_upset(sites: pd.DataFrame, field: MembershipField, output_path, label: str, pass_only: bool = False) -> None:
    fig = plt.figure(figsize=double_column_figsize(3.6))
    series = _build_upset_series(sites, field, pass_only=pass_only)
    if series is None or series.empty:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No variants", ha="center", va="center")
        ax.set_axis_off()
        save_figure(fig, output_path)
        return
    upset = CompatibleUpSet(
        series,
        orientation="horizontal",
        sort_by="cardinality",
        sort_categories_by=None,
        min_subset_size=None,
        max_subset_size=None,
        max_subset_rank=None,
        include_empty_subsets=False,
    )
    upset.plot(fig=fig)
    fig.suptitle(f"{field.title} observed combinations: {label}", y=0.98)
    save_figure(fig, output_path)


class UpSetModule(AnalysisModule):
    @property
    def name(self) -> str:
        return "upset"

    def run(self, data: AggregatedData, config: AnalysisConfig) -> None:
        output_dir = self.output_dir(config)
        tables_dir = output_dir / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)
        for label, sites in ((data.label_a, data.sites_a), (data.label_b, data.sites_b)):
            for field in _MEMBERSHIP_FIELDS:
                table = build_membership_count_table(sites, field.name, pass_only=config.pass_only)
                write_tsv_gz(table, tables_dir / f"{field.name}_combinations.{label}.tsv")
                _plot_upset(sites, field, output_dir / f"{field.name}.upset.{label}.png", label, pass_only=config.pass_only)
