import os
import re
import json
from pathlib import Path
from IPython.display import display, Image, clear_output
import ipywidgets as widgets
from datetime import datetime

# Subtypes for "correct" classifications
CORRECT_SUBTYPES = [
    "typical",
    "atypical",
    "homozygous",
    "complex variant",
    "mosaic",
    "other"
]


# Display review statistics
def display_statistics(manifest_file):
    """Display statistics from the manifest."""
    if not os.path.exists(manifest_file):
        print("No manifest file found")
        return
    
    with open(manifest_file, 'r') as f:
        manifest = json.load(f)
    
    total_reviews = len(manifest)
    correct_count = sum(1 for r in manifest.values() if r['classification'] == 'correct')
    incorrect_count = total_reviews - correct_count
    
    print(f"Total reviews: {total_reviews}")
    print(f"Correct: {correct_count} ({correct_count/total_reviews*100:.1f}%)")
    print(f"Incorrect: {incorrect_count} ({incorrect_count/total_reviews*100:.1f}%)")
    print()
    
    # Subtype breakdown
    if correct_count > 0:
        print("Correct subtypes:")
        subtype_counts = {}
        for review in manifest.values():
            if review['classification'] == 'correct' and 'subtype' in review:
                subtype = review['subtype']
                subtype_counts[subtype] = subtype_counts.get(subtype, 0) + 1
        
        for subtype, count in sorted(subtype_counts.items(), key=lambda x: -x[1]):
            print(f"  {subtype}: {count} ({count/correct_count*100:.1f}%)")


def parse_filename(filename):
    """
    Parse plot filename to extract components.
    Format: chr_start_end_{GD_ID}_{svtype}_{SAMPLE_ID}.jpg
    Also handles duplicated coordinates in filename: chr_start_end_{GD_ID}_{svtype}_chr_start_end_{svtype}_{SAMPLE_ID}
    """
    # Remove extension
    base = filename.replace('.jpg', '').replace('.jpeg', '').replace('.png', '')
    
    # Try to parse using regex - handle both underscore and dash as separators
    # New Pattern: chr{N}_{start}_{end}_{GD_ID}_{svtype}_{SAMPLE_ID}
    # Also handles: chr{N}_{start}-{end}_{GD_ID}_...
    pattern = r'(chr[^_]+)_(\d+)[-_](\d+)_(GD_[^_]+(?:_[^_]+)*?)_(DEL|DUP)_(.+)'
    match = re.match(pattern, base)
    
    if match:
        info = {
            'chr': match.group(1),
            'start': int(match.group(2)),
            'end': int(match.group(3)),
            'gd_id': match.group(4),
            'svtype': match.group(5),
            'sample_id': match.group(6),
            'filename': filename
        }
        
        # Check for redundant coordinate/svtype prefix in sample_id
        # Example filename: chr1_100_200_GD_ID_DEL_chr1_100_200_DEL_sample_id.jpg
        # Extracted sample_id would be: chr1_100_200_DEL_sample_id
        
        # Construct expected prefix
        prefix = f"{info['chr']}_{info['start']}_{info['end']}_{info['svtype']}_"
        if info['sample_id'].startswith(prefix):
            info['sample_id'] = info['sample_id'][len(prefix):]
            
        return info
    else:
        return None

def load_manifest(manifest_file):
    """Load existing manifest or create new one."""
    if os.path.exists(manifest_file):
        with open(manifest_file, 'r') as f:
            return json.load(f)
    return {}

def save_manifest(manifest, manifest_file):
    """Save manifest to file."""
    with open(manifest_file, 'w') as f:
        json.dump(manifest, f, indent=2)

def get_plot_files(directory):
    """Get all plot files from directory."""
    plot_dir = Path(directory)
    if not plot_dir.exists():
        print(f"Warning: Directory {directory} does not exist")
        return []
    
    extensions = ['.jpg', '.jpeg', '.png']
    files = []
    for ext in extensions:
        files.extend(plot_dir.glob(f'*{ext}'))
    
    return sorted([str(f) for f in files])

class PlotReviewApp:
    def __init__(self, plot_directory, manifest_file, gd_regions_df=None):
        self.plot_directory = plot_directory
        self.manifest_file = manifest_file
        self.manifest = load_manifest(self.manifest_file)

        if gd_regions_df is not None:
            self.gd_regions_df = gd_regions_df
        else:
            # Fallback if variable is available globally, though cleaner to pass it in
            self.gd_regions_df = None

        self.gd_meta_by_key = self._build_gd_meta_lookup(self.gd_regions_df)
        self.all_plot_files = get_plot_files(plot_directory)
        self.plot_files = self.all_plot_files.copy()
        self.current_index = 0

        # Find first unreviewed plot
        for i, plot_file in enumerate(self.plot_files):
            if plot_file not in self.manifest:
                self.current_index = i
                break

        self.setup_ui()

    @staticmethod
    def _normalize_chr(chrom):
        if chrom is None:
            return None
        chrom = str(chrom)
        if chrom.startswith("chr"):
            return chrom
        return f"chr{chrom}"

    @staticmethod
    def _find_first_column(df, candidates):
        if df is None:
            return None
        for col in candidates:
            if col in df.columns:
                return col
        # fallback: case-insensitive
        lower_to_original = {c.lower(): c for c in df.columns}
        for col in candidates:
            if col.lower() in lower_to_original:
                return lower_to_original[col.lower()]
        return None

    @staticmethod
    def _unique_str_options(series):
        if series is None:
            return ["All"]
        values = []
        for v in series.dropna().unique().tolist():
            s = str(v).strip()
            if s:
                values.append(s)
        return ["All"] + sorted(set(values))

    def _build_gd_meta_lookup(self, gd_regions_df):
        """Build lookup from GD table keyed by plot_info fields (chr, start, end, svtype)."""
        if gd_regions_df is None or not hasattr(gd_regions_df, "columns"):
            return {}

        chr_col = self._find_first_column(gd_regions_df, ["chr", "chrom", "chromosome"])
        start_col = self._find_first_column(gd_regions_df, ["start_GRCh38", "start", "start_hg38"])
        end_col = self._find_first_column(gd_regions_df, ["end_GRCh38", "end", "end_hg38"])
        # gd_id_col = self._find_first_column(gd_regions_df, ["GD_ID", "gd_id", "gdid"]) # Not using GD_ID for lookup key due to format mismatch
        svtype_col = self._find_first_column(gd_regions_df, ["svtype", "SVTYPE"])

        cluster_col = self._find_first_column(gd_regions_df, ["cluster", "cluster_ID", "cluster_id"])
        terminal_col = self._find_first_column(gd_regions_df, ["terminal", "Terminal", "is_terminal"])
        nahr_col = self._find_first_column(gd_regions_df, ["NAHR", "nahr"])

        required = [chr_col, start_col, end_col, svtype_col]
        if any(c is None for c in required):
            return {}

        lookup = {}
        for _, row in gd_regions_df.iterrows():
            try:
                chrom = self._normalize_chr(row[chr_col])
                start = int(row[start_col])
                end = int(row[end_col])
                # gd_id = str(row[gd_id_col])
                svtype = str(row[svtype_col])
            except Exception:
                continue

            meta = {
                "cluster": None if cluster_col is None else row.get(cluster_col),
                "terminal": None if terminal_col is None else row.get(terminal_col),
                "NAHR": None if nahr_col is None else row.get(nahr_col),
                "svtype": svtype,
            }

            # key = (chrom, start, end, gd_id, svtype)
            # Use reduced key to avoid GD_ID string mismatches
            key = (chrom, start, end, svtype)
            lookup[key] = meta
        return lookup

    def _make_plot_key(self, plot_info):
        if not plot_info:
            return None
        return (
            self._normalize_chr(plot_info.get("chr")),
            int(plot_info.get("start")),
            int(plot_info.get("end")),
            # str(plot_info.get("gd_id")), # inconsistent format
            str(plot_info.get("svtype")),
        )

    def get_unique_gd_ids(self):
        """Extract unique GD IDs from all plot files."""
        gd_ids = set()
        for plot_file in self.all_plot_files:
            plot_info = parse_filename(os.path.basename(plot_file))
            if plot_info and "gd_id" in plot_info:
                gd_ids.add(plot_info["gd_id"])
        return ["All"] + sorted(list(gd_ids))

    def apply_filters(self, change=None):
        """Apply filters to plot files."""
        gd_filter = self.gd_filter_dropdown.value
        sample_filter = self.sample_filter_text.value.strip()
        review_status = self.review_status_dropdown.value
        classification_filter = self.classification_filter_dropdown.value

        cluster_filter = self.cluster_filter_dropdown.value
        terminal_filter = self.terminal_filter_dropdown.value
        nahr_filter = self.nahr_filter_dropdown.value
        svtype_filter = self.svtype_filter_dropdown.value

        # Start with all files
        filtered_files = self.all_plot_files.copy()

        # Apply GD ID filter
        if gd_filter != "All":
            temp_files = []
            for plot_file in filtered_files:
                plot_info = parse_filename(os.path.basename(plot_file))
                if plot_info and plot_info.get("gd_id") == gd_filter:
                    temp_files.append(plot_file)
            filtered_files = temp_files

        # Apply Sample ID filter (case-insensitive substring match)
        if sample_filter:
            temp_files = []
            for plot_file in filtered_files:
                plot_info = parse_filename(os.path.basename(plot_file))
                if plot_info and sample_filter.lower() in plot_info.get("sample_id", "").lower():
                    temp_files.append(plot_file)
            filtered_files = temp_files

        # Apply GD-table region filters (cluster / terminal / NAHR / svtype)
        if (
            cluster_filter != "All"
            or terminal_filter != "All"
            or nahr_filter != "All"
            or svtype_filter != "All"
        ):
            temp_files = []
            for plot_file in filtered_files:
                plot_info = parse_filename(os.path.basename(plot_file))
                if not plot_info:
                    continue

                if svtype_filter != "All" and plot_info.get("svtype") != svtype_filter:
                    continue

                plot_key = self._make_plot_key(plot_info)
                meta = self.gd_meta_by_key.get(plot_key)

                if meta is None:
                    # If any of the metadata-dependent filters are active, require a GD table match
                    if cluster_filter != "All" or terminal_filter != "All" or nahr_filter != "All":
                        continue
                    temp_files.append(plot_file)
                    continue

                if cluster_filter != "All" and str(meta.get("cluster")).strip() != cluster_filter:
                    continue
                if terminal_filter != "All" and str(meta.get("terminal")).strip() != terminal_filter:
                    continue
                if nahr_filter != "All" and str(meta.get("NAHR")).strip() != nahr_filter:
                    continue

                temp_files.append(plot_file)
            filtered_files = temp_files

        # Apply Review Status filter
        if review_status != "All":
            temp_files = []
            for plot_file in filtered_files:
                is_reviewed = plot_file in self.manifest
                if review_status == "Reviewed" and is_reviewed:
                    temp_files.append(plot_file)
                elif review_status == "Not Reviewed" and not is_reviewed:
                    temp_files.append(plot_file)
            filtered_files = temp_files

        # Apply Classification filter (only meaningful for reviewed plots)
        if classification_filter != "All":
            temp_files = []
            for plot_file in filtered_files:
                if plot_file not in self.manifest:
                    continue
                classification = self.manifest[plot_file].get("classification")
                if classification_filter == "Correct" and classification == "correct":
                    temp_files.append(plot_file)
                elif classification_filter == "Incorrect" and classification == "incorrect":
                    temp_files.append(plot_file)
            filtered_files = temp_files

        self.plot_files = filtered_files

        # Reset to first file in filtered list
        if self.plot_files:
            self.current_index = 0

        self.update_display()

    def clear_filters(self, b=None):
        """Clear all filters."""
        self.gd_filter_dropdown.value = "All"
        self.sample_filter_text.value = ""
        self.review_status_dropdown.value = "All"
        self.classification_filter_dropdown.value = "All"

        self.cluster_filter_dropdown.value = "All"
        self.terminal_filter_dropdown.value = "All"
        self.nahr_filter_dropdown.value = "All"
        self.svtype_filter_dropdown.value = "All"

        self.apply_filters()

    def setup_ui(self):
        """Setup the user interface."""
        # Header
        self.header = widgets.HTML()

        # Filters
        self.gd_filter_dropdown = widgets.Dropdown(
            options=self.get_unique_gd_ids(),
            value="All",
            description="Filter GD:",
            layout=widgets.Layout(width="400px"),
        )

        self.sample_filter_text = widgets.Text(
            placeholder="Enter sample ID to filter...",
            description="Filter Sample:",
            layout=widgets.Layout(width="400px"),
            continuous_update=False
        )

        self.review_status_dropdown = widgets.Dropdown(
            options=["All", "Reviewed", "Not Reviewed"],
            value="All",
            description="Review Status:",
            layout=widgets.Layout(width="400px"),
        )

        self.classification_filter_dropdown = widgets.Dropdown(
            options=["All", "Correct", "Incorrect"],
            value="All",
            description="Classification:",
            layout=widgets.Layout(width="400px"),
        )

        # GD table-driven filters
        df = self.gd_regions_df if hasattr(self.gd_regions_df, "columns") else None
        cluster_col = self._find_first_column(df, ["cluster", "cluster_ID", "cluster_id"]) if df is not None else None
        terminal_col = self._find_first_column(df, ["terminal", "Terminal", "is_terminal"]) if df is not None else None
        nahr_col = self._find_first_column(df, ["NAHR", "nahr"]) if df is not None else None
        svtype_col = self._find_first_column(df, ["svtype", "SVTYPE"]) if df is not None else None

        cluster_options = self._unique_str_options(df[cluster_col]) if (df is not None and cluster_col is not None) else ["All"]
        terminal_options = self._unique_str_options(df[terminal_col]) if (df is not None and terminal_col is not None) else ["All"]
        nahr_options = self._unique_str_options(df[nahr_col]) if (df is not None and nahr_col is not None) else ["All"]
        svtype_options = self._unique_str_options(df[svtype_col]) if (df is not None and svtype_col is not None) else ["All", "DEL", "DUP"]

        self.cluster_filter_dropdown = widgets.Dropdown(
            options=cluster_options,
            value="All",
            description="Cluster:",
            layout=widgets.Layout(width="400px"),
            disabled=(len(cluster_options) <= 1),
        )

        self.terminal_filter_dropdown = widgets.Dropdown(
            options=terminal_options,
            value="All",
            description="Terminal:",
            layout=widgets.Layout(width="400px"),
            disabled=(len(terminal_options) <= 1),
        )

        self.nahr_filter_dropdown = widgets.Dropdown(
            options=nahr_options,
            value="All",
            description="NAHR:",
            layout=widgets.Layout(width="400px"),
            disabled=(len(nahr_options) <= 1),
        )

        self.svtype_filter_dropdown = widgets.Dropdown(
            options=svtype_options,
            value="All",
            description="SVTYPE:",
            layout=widgets.Layout(width="400px"),
        )

        self.apply_filter_button = widgets.Button(
            description="Apply Filters",
            button_style="warning",
            icon="filter",
        )

        self.clear_filter_button = widgets.Button(
            description="Clear Filters",
            button_style="",
            icon="times",
        )

        # Image display
        self.image_widget = widgets.Output()

        # Plot info
        self.info_widget = widgets.HTML()

        # Quick review buttons
        self.quick_correct_typical_button = widgets.Button(
            description="Correct typical",
            button_style="success",
            icon="check",
            layout=widgets.Layout(width='auto', min_width='150px')
        )
        self.quick_incorrect_button = widgets.Button(
            description="Incorrect",
            button_style="danger",
            icon="times",
            layout=widgets.Layout(width='auto', min_width='150px')
        )
        self.correct_in_cluster_button = widgets.Button(
            description="Pick correct for sample cluster",
            button_style="primary",
            icon="check-circle",
            layout=widgets.Layout(width='auto', min_width='230px')
        )
        self.filter_to_sample_button = widgets.Button(
            description="Filter to sample",
            button_style="warning",
            icon="filter",
            layout=widgets.Layout(width='auto', min_width='200px')
        )
        self.filter_to_region_button = widgets.Button(
            description="Filter to region",
            button_style="warning",
            icon="filter",
            layout=widgets.Layout(width='auto', min_width='200px')
        )
        self.filter_to_cluster_button = widgets.Button(
            description="Filter to cluster",
            button_style="warning",
            icon="filter",
            layout=widgets.Layout(width='auto', min_width='200px')
        )

        # Classification controls
        self.classification_radio = widgets.RadioButtons(
            options=["correct", "incorrect"],
            description="Classification:",
            disabled=False,
        )

        self.subtype_dropdown = widgets.Dropdown(
            options=CORRECT_SUBTYPES,
            description="Subtype:",
            disabled=False,
        )

        # Notes field
        self.notes_text = widgets.Textarea(
            placeholder="Optional notes...",
            description="Notes:",
            layout=widgets.Layout(width="600px", height="80px"),
        )

        # Navigation buttons
        self.prev_button = widgets.Button(
            description="← Previous",
            button_style="info",
            disabled=True,
        )
        self.next_button = widgets.Button(
            description="Next →",
            button_style="info",
        )
        self.submit_button = widgets.Button(
            description="Submit Review",
            button_style="success",
            icon="check",
        )
        self.clear_review_button = widgets.Button(
            description="Clear Review",
            button_style="danger",
            icon="trash",
        )

        # Progress info
        self.progress_widget = widgets.HTML()

        # Event handlers
        self.prev_button.on_click(self.on_prev)
        self.next_button.on_click(self.on_next)
        self.submit_button.on_click(self.on_submit)
        self.clear_review_button.on_click(self.on_clear_review)
        self.quick_correct_typical_button.on_click(self.on_quick_correct_typical)
        self.quick_incorrect_button.on_click(self.on_quick_incorrect)
        self.filter_to_sample_button.on_click(self.on_filter_to_sample)
        self.filter_to_region_button.on_click(self.on_filter_to_region)
        self.filter_to_cluster_button.on_click(self.on_filter_to_cluster)
        self.correct_in_cluster_button.on_click(self.on_correct_in_cluster)
        self.classification_radio.observe(self.on_classification_change, names="value")
        self.apply_filter_button.on_click(self.apply_filters)
        self.clear_filter_button.on_click(self.clear_filters)

        # Layout
        quick_box = widgets.VBox(
            [
                widgets.HTML("<b>Quick review:</b>"),
                widgets.HBox([self.quick_correct_typical_button, self.quick_incorrect_button, self.correct_in_cluster_button]),
                widgets.HBox([self.filter_to_sample_button, self.filter_to_region_button, self.filter_to_cluster_button]),
            ]
        )

        nav_buttons = widgets.HBox([self.prev_button, self.next_button])
        action_buttons = widgets.HBox([self.submit_button, self.clear_review_button])
        
        # Review controls
        self.review_controls = widgets.VBox(
            [
                nav_buttons,
                quick_box,
                self.classification_radio,
                self.subtype_dropdown,
                self.notes_text,
                action_buttons,
            ]
        )
        
        # Filters
        self.filter_box = widgets.VBox(
            [
                widgets.HTML("<h3>Filters</h3>"),
                self.gd_filter_dropdown,
                self.sample_filter_text,
                widgets.HTML("<b>Region (GD table)</b>"),
                self.cluster_filter_dropdown,
                self.terminal_filter_dropdown,
                self.nahr_filter_dropdown,
                self.svtype_filter_dropdown,
                self.review_status_dropdown,
                self.classification_filter_dropdown,
                widgets.HBox([self.apply_filter_button, self.clear_filter_button]),
            ]
        )
        
        # Combine review controls and filters vertically
        self.controls_box = widgets.VBox(
            [
                self.review_controls,
                widgets.HTML("<hr>"),
                self.filter_box,
            ],
            layout=widgets.Layout(width='100%')
        )

        self.container = widgets.VBox(
            [
                self.header,
                widgets.HTML("<hr>"),
                self.progress_widget,
                self.info_widget,
                self.image_widget,
                widgets.HTML("<hr>"),
                self.controls_box,
            ]
        )

    def on_quick_correct_typical(self, b):
        """One-click submit: correct + typical."""
        if not self.plot_files:
            return
        self.classification_radio.value = "correct"
        if "typical" in CORRECT_SUBTYPES:
            self.subtype_dropdown.value = "typical"
        else:
            self.subtype_dropdown.value = CORRECT_SUBTYPES[0]
        self.on_submit(b)

    def on_quick_incorrect(self, b):
        """One-click submit: incorrect."""
        if not self.plot_files:
            return
        self.classification_radio.value = "incorrect"
        self.on_submit(b)

    def on_filter_to_sample(self, b):
        """Filter to current sample."""
        if not self.plot_files:
            return
        
        current_file = self.plot_files[self.current_index]
        plot_info = parse_filename(os.path.basename(current_file))
        
        if plot_info and "sample_id" in plot_info:
            self.sample_filter_text.value = plot_info["sample_id"]
            self.apply_filters()

    def on_filter_to_region(self, b):
        """Filter to current region (GD_ID)."""
        if not self.plot_files:
            return
        
        current_file = self.plot_files[self.current_index]
        plot_info = parse_filename(os.path.basename(current_file))
        
        if plot_info and "gd_id" in plot_info:
            self.gd_filter_dropdown.value = plot_info["gd_id"]
            self.apply_filters()

    def on_filter_to_cluster(self, b):
        """Filter to current cluster."""
        if not self.plot_files:
            return
        
        current_file = self.plot_files[self.current_index]
        plot_info = parse_filename(os.path.basename(current_file))
        
        if not plot_info:
            return
        
        plot_key = self._make_plot_key(plot_info)
        meta = self.gd_meta_by_key.get(plot_key)
        
        if meta and meta.get("cluster") is not None:
            cluster_value = str(meta.get("cluster")).strip()
            if cluster_value in self.cluster_filter_dropdown.options:
                self.cluster_filter_dropdown.value = cluster_value
                self.apply_filters()

    def _apply_filters_with_override(self, override_filter=None, override_value=None):
        """Apply current filters with one filter overridden. Returns count of matching files."""
        # Get current filter values
        gd_filter = self.gd_filter_dropdown.value
        sample_filter = self.sample_filter_text.value.strip()
        review_status = self.review_status_dropdown.value
        classification_filter = self.classification_filter_dropdown.value
        cluster_filter = self.cluster_filter_dropdown.value
        terminal_filter = self.terminal_filter_dropdown.value
        nahr_filter = self.nahr_filter_dropdown.value
        svtype_filter = self.svtype_filter_dropdown.value
        
        # Override one filter if specified
        if override_filter == 'sample':
            sample_filter = override_value
        elif override_filter == 'region':
            gd_filter = override_value
        elif override_filter == 'cluster':
            cluster_filter = override_value
        
        # Start with all files
        filtered_files = self.all_plot_files.copy()

        # Apply GD ID filter
        if gd_filter != "All":
            temp_files = []
            for plot_file in filtered_files:
                plot_info = parse_filename(os.path.basename(plot_file))
                if plot_info and plot_info.get("gd_id") == gd_filter:
                    temp_files.append(plot_file)
            filtered_files = temp_files

        # Apply Sample ID filter (case-insensitive substring match)
        if sample_filter:
            temp_files = []
            for plot_file in filtered_files:
                plot_info = parse_filename(os.path.basename(plot_file))
                if plot_info and sample_filter.lower() in plot_info.get("sample_id", "").lower():
                    temp_files.append(plot_file)
            filtered_files = temp_files

        # Apply GD-table region filters (cluster / terminal / NAHR / svtype)
        if (
            cluster_filter != "All"
            or terminal_filter != "All"
            or nahr_filter != "All"
            or svtype_filter != "All"
        ):
            temp_files = []
            for plot_file in filtered_files:
                plot_info = parse_filename(os.path.basename(plot_file))
                if not plot_info:
                    continue

                if svtype_filter != "All" and plot_info.get("svtype") != svtype_filter:
                    continue

                plot_key = self._make_plot_key(plot_info)
                meta = self.gd_meta_by_key.get(plot_key)

                if meta is None:
                    if cluster_filter != "All" or terminal_filter != "All" or nahr_filter != "All":
                        continue
                    temp_files.append(plot_file)
                    continue

                if cluster_filter != "All" and str(meta.get("cluster")).strip() != cluster_filter:
                    continue
                if terminal_filter != "All" and str(meta.get("terminal")).strip() != terminal_filter:
                    continue
                if nahr_filter != "All" and str(meta.get("NAHR")).strip() != nahr_filter:
                    continue

                temp_files.append(plot_file)
            filtered_files = temp_files

        # Apply Review Status filter
        if review_status != "All":
            temp_files = []
            for plot_file in filtered_files:
                is_reviewed = plot_file in self.manifest
                if review_status == "Reviewed" and is_reviewed:
                    temp_files.append(plot_file)
                elif review_status == "Not Reviewed" and not is_reviewed:
                    temp_files.append(plot_file)
            filtered_files = temp_files

        # Apply Classification filter
        if classification_filter != "All":
            temp_files = []
            for plot_file in filtered_files:
                if plot_file not in self.manifest:
                    continue
                classification = self.manifest[plot_file].get("classification")
                if classification_filter == "Correct" and classification == "correct":
                    temp_files.append(plot_file)
                elif classification_filter == "Incorrect" and classification == "incorrect":
                    temp_files.append(plot_file)
            filtered_files = temp_files
        
        return len(filtered_files)
    
    def _calculate_filter_counts(self, plot_info):
        """Calculate the number of plots that would result from each filter type."""
        counts = {
            'sample': 0,
            'region': 0,
            'cluster': 0
        }
        
        if not plot_info:
            return counts
        
        sample_id = plot_info.get('sample_id')
        gd_id = plot_info.get('gd_id')
        
        # Count plots for sample filter
        if sample_id:
            counts['sample'] = self._apply_filters_with_override('sample', sample_id)
        
        # Count plots for region filter
        if gd_id:
            counts['region'] = self._apply_filters_with_override('region', gd_id)
        
        # Count plots for cluster filter
        plot_key = self._make_plot_key(plot_info)
        meta = self.gd_meta_by_key.get(plot_key)
        if meta and meta.get('cluster') is not None:
            cluster_value = str(meta.get('cluster')).strip()
            counts['cluster'] = self._apply_filters_with_override('cluster', cluster_value)
        
        return counts
    
    def _get_other_plots_in_cluster(self, current_file, plot_info):
        """Get other plots for the same sample and cluster."""
        if not plot_info:
            return []
        
        sample_id = plot_info.get("sample_id")
        plot_key = self._make_plot_key(plot_info)
        meta = self.gd_meta_by_key.get(plot_key)
        
        if not sample_id or not meta:
            return []
        
        cluster_id = meta.get("cluster")
        if cluster_id is None:
            return []
        
        # Find all plots with same sample_id and cluster
        other_plots = []
        for plot_file in self.all_plot_files:
            if plot_file == current_file:
                continue
            
            file_info = parse_filename(os.path.basename(plot_file))
            if not file_info:
                continue
            
            if file_info.get("sample_id") != sample_id:
                continue
            
            file_key = self._make_plot_key(file_info)
            file_meta = self.gd_meta_by_key.get(file_key)
            
            if file_meta and str(file_meta.get("cluster")) == str(cluster_id):
                other_plots.append(plot_file)
        
        return other_plots

    def on_correct_in_cluster(self, b):
        """Mark current plot as correct and all others in sample+cluster as incorrect."""
        if not self.plot_files:
            return
        
        current_file = self.plot_files[self.current_index]
        plot_info = parse_filename(os.path.basename(current_file))
        
        if not plot_info:
            return
        
        # Get other plots in the same sample+cluster
        other_plots = self._get_other_plots_in_cluster(current_file, plot_info)
        
        if len(other_plots) == 0:
            return
        
        # Mark all other plots as incorrect
        for plot_file in other_plots:
            other_plot_info = parse_filename(os.path.basename(plot_file))
            
            review = {
                "classification": "incorrect",
                "notes": "Auto-marked as incorrect (another plot in cluster marked as correct)",
                "timestamp": datetime.now().isoformat(),
                "filename": plot_file,
            }
            
            if other_plot_info:
                review.update(other_plot_info)
            
            self.manifest[plot_file] = review
        
        # Save manifest
        save_manifest(self.manifest, self.manifest_file)
        
        # Mark current plot as correct typical
        self.classification_radio.value = "correct"
        if "typical" in CORRECT_SUBTYPES:
            self.subtype_dropdown.value = "typical"
        else:
            self.subtype_dropdown.value = CORRECT_SUBTYPES[0]
        self.on_submit(b)

    def on_classification_change(self, change):
        """Handle classification radio button change."""
        if change["new"] == "correct":
            self.subtype_dropdown.disabled = False
        else:
            self.subtype_dropdown.disabled = True

    def on_clear_review(self, b):
        """Clear review for current plot."""
        if not self.plot_files:
            return

        current_file = self.plot_files[self.current_index]

        if current_file in self.manifest:
            del self.manifest[current_file]
            save_manifest(self.manifest, self.manifest_file)
            self.update_display()

    def update_display(self):
        """Update the display with current plot."""
        if not self.plot_files:
            self.header.value = "<h2>Genomic Disorder Plot Review</h2>"
            self.progress_widget.value = ""
            self.info_widget.value = ""
            with self.image_widget:
                clear_output(wait=True)
                display(widgets.HTML("<h2 style='color: orange;'>No plots match current filters!</h2>"))
            # Hide review controls but keep filters visible
            self.review_controls.layout.display = 'none'
            return
        
        # Show all controls
        self.review_controls.layout.display = 'block'

        if self.current_index >= len(self.plot_files):
            self.header.value = "<h2 style='color: green;'>All filtered plots reviewed!</h2>"
            return

        current_file = self.plot_files[self.current_index]

        # Update header
        self.header.value = "<h2>Genomic Disorder Plot Review</h2>"

        # Update progress (based on ALL files, not filtered)
        total_reviewed_count = sum(1 for f in self.all_plot_files if f in self.manifest)
        filtered_reviewed_count = sum(1 for f in self.plot_files if f in self.manifest)

        filter_status = ""
        if len(self.plot_files) < len(self.all_plot_files):
            filter_status = (
                f" | <b>Filtered:</b> {len(self.plot_files)} of {len(self.all_plot_files)} plots "
                f"({filtered_reviewed_count} reviewed)"
            )

        self.progress_widget.value = f"""
        <div style='background-color: #e8f4f8; padding: 10px; border-radius: 5px;'>
            <b>Overall Progress:</b> {total_reviewed_count}/{len(self.all_plot_files)} reviewed
            ({total_reviewed_count/len(self.all_plot_files)*100:.1f}%){filter_status} |
            <b>Current:</b> {self.current_index + 1}/{len(self.plot_files)}
        </div>
        """

        # Determine review status string
        if current_file in self.manifest:
            review = self.manifest[current_file]
            classification = review["classification"]
            subtype_str = ""
            if classification == "correct":
                subtype_str = f" ({review.get('subtype', 'N/A')})"
            
            # Format timestamp
            timestamp_str = ""
            if "timestamp" in review:
                try:
                    # Parse ISO format timestamp and format it nicely
                    from datetime import datetime as dt
                    ts = dt.fromisoformat(review["timestamp"])
                    timestamp_str = f" at {ts.strftime('%Y-%m-%d %H:%M:%S')}"
                except:
                    timestamp_str = f" at {review['timestamp']}"
            
            status_html = f"<br><b>Status:</b> <span style='color: green;'>Reviewed</span> - {classification.upper()}{subtype_str}{timestamp_str}"
        else:
            status_html = "<br><b>Status:</b> <span style='color: red;'>Not Reviewed</span>"

        # Parse filename
        plot_info = parse_filename(os.path.basename(current_file))
        if plot_info:
            self.info_widget.value = f"""
            <div style='background-color: #f0f0f0; padding: 10px; border-radius: 5px; margin: 10px 0;'>
                <b>Sample ID:</b> {plot_info['sample_id']}<br>
                <b>GD ID:</b> {plot_info['gd_id']}<br>
                <b>Region:</b> {plot_info['chr']}:{plot_info['start']}-{plot_info['end']}<br>
                <b>SV Type:</b> {plot_info['svtype']}<br>
                <b>File:</b> {os.path.basename(current_file)}
                {status_html}
            </div>
            """
        else:
            self.info_widget.value = f"""
            <div style='background-color: #fff3cd; padding: 10px; border-radius: 5px; margin: 10px 0;'>
                <b>Warning:</b> Could not parse filename<br>
                <b>File:</b> {os.path.basename(current_file)}
                {status_html}
            </div>
            """

        # Display image
        with self.image_widget:
            clear_output(wait=True)
            if os.path.exists(current_file):
                display(Image(filename=current_file, width=900))
            else:
                print(f"Image not found: {current_file}")

        # Load existing review if present
        if current_file in self.manifest:
            review = self.manifest[current_file]
            self.classification_radio.value = review["classification"]
            if review["classification"] == "correct":
                self.subtype_dropdown.value = review.get("subtype", CORRECT_SUBTYPES[0])
                self.subtype_dropdown.disabled = False
            else:
                self.subtype_dropdown.disabled = True
            self.notes_text.value = review.get("notes", "")
        else:
            # Reset to defaults
            self.classification_radio.value = "correct"
            self.subtype_dropdown.value = CORRECT_SUBTYPES[0]
            self.subtype_dropdown.disabled = False
            self.notes_text.value = ""

        # Update button states
        self.prev_button.disabled = self.current_index == 0
        self.next_button.disabled = self.current_index >= len(self.plot_files) - 1
        
        # Update correct_in_cluster_button state
        if plot_info:
            other_cluster_plots = self._get_other_plots_in_cluster(current_file, plot_info)
            self.correct_in_cluster_button.disabled = len(other_cluster_plots) == 0
        else:
            self.correct_in_cluster_button.disabled = True
        
        # Update filter button descriptions with counts
        if plot_info:
            filter_counts = self._calculate_filter_counts(plot_info)
            self.filter_to_sample_button.description = f"Filter to sample ({filter_counts['sample']})"
            self.filter_to_region_button.description = f"Filter to region ({filter_counts['region']})"
            self.filter_to_cluster_button.description = f"Filter to cluster ({filter_counts['cluster']})"
        else:
            self.filter_to_sample_button.description = "Filter to sample"
            self.filter_to_region_button.description = "Filter to region"
            self.filter_to_cluster_button.description = "Filter to cluster"

    def on_prev(self, b):
        """Navigate to previous plot."""
        if self.current_index > 0:
            self.current_index -= 1
            self.update_display()

    def on_next(self, b):
        """Navigate to next plot."""
        if self.current_index < len(self.plot_files) - 1:
            self.current_index += 1
            self.update_display()

    def on_submit(self, b):
        """Submit current review."""
        if not self.plot_files:
            return

        current_file = self.plot_files[self.current_index]
        plot_info = parse_filename(os.path.basename(current_file))

        # Create review record
        review = {
            "classification": self.classification_radio.value,
            "notes": self.notes_text.value,
            "timestamp": datetime.now().isoformat(),
            "filename": current_file,
        }

        if self.classification_radio.value == "correct":
            review["subtype"] = self.subtype_dropdown.value

        if plot_info:
            review.update(plot_info)

        # Save to manifest
        self.manifest[current_file] = review
        save_manifest(self.manifest, self.manifest_file)

        # Move to next unreviewed plot
        found_next = False
        for i in range(self.current_index + 1, len(self.plot_files)):
            if self.plot_files[i] not in self.manifest:
                self.current_index = i
                found_next = True
                break

        if not found_next and self.current_index < len(self.plot_files) - 1:
            self.current_index += 1

        self.update_display()

    def display(self):
        """Display the app."""
        self.update_display()
        display(self.container)
