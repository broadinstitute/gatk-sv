#!/bin/python

import argparse
import json
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

from sklearn.metrics import confusion_matrix, classification_report


def read_file_list(file_list_path):
    """
    Read a list of file paths from a text file.

    Args:
        file_list_path: Path to text file with file paths (one per line)

    Returns:
        List of file paths as strings
    """
    with open(file_list_path, "r") as f:
        file_paths = [line.strip() for line in f if line.strip()]
    return file_paths


def concatenate_tsvs(file_paths):
    """
    Read and concatenate multiple TSV files into a single DataFrame.

    Args:
        file_paths: List of paths to TSV files

    Returns:
        Concatenated pandas DataFrame
    """
    dfs = []
    print(f"Reading {len(file_paths)} TSV files...")
    for file_path in file_paths:
        if os.path.exists(file_path):
            df = pd.read_csv(file_path, sep="\t")
            dfs.append(df)
        else:
            print(f"Warning: File not found: {file_path}")

    if not dfs:
        raise ValueError("No valid TSV files found to concatenate")

    combined_df = pd.concat(dfs, ignore_index=True)
    print(f"Combined data: {len(combined_df)} rows, {len(combined_df.columns)} columns")
    return combined_df


def get_chromosome_type(chrom):
    """
    Classify chromosome as autosomal, chrX, or chrY.

    Args:
        chrom: Chromosome name (e.g., 'chr1', 'chrX', 'chrY')

    Returns:
        String: 'chrX', 'chrY', or 'Autosomal'
    """
    if chrom == "chrX":
        return "chrX"
    elif chrom == "chrY":
        return "chrY"
    else:
        return "Autosomal"


def format_column_name(col):
    """
    Format column name for display.

    Args:
        col: Column name string

    Returns:
        Formatted column name with underscores replaced by spaces and title cased
    """
    return col.replace("_", " ").title()


def save_and_close_plot(output_dir, filename):
    """
    Save the current matplotlib plot and close it.

    Args:
        output_dir: Directory path where plot will be saved
        filename: Name of output file

    Returns:
        None
    """
    plt.tight_layout()
    diagnostics_dir = os.path.join(output_dir, "diagnostics")
    os.makedirs(diagnostics_dir, exist_ok=True)
    plt.savefig(os.path.join(diagnostics_dir, filename), dpi=300)
    plt.close()
    print(f"  Created histogram: diagnostics/{filename}")


def plot_single_histogram(data, x_col, title, xlabel, output_dir, filename, bins=50):
    """
    Create and save a single histogram plot.

    Args:
        data: Pandas DataFrame containing the data
        x_col: Column name to plot on x-axis
        title: Plot title
        xlabel: X-axis label
        output_dir: Directory path where plot will be saved
        filename: Name of output file
        bins: Number of histogram bins (default: 50)

    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(
        data=data,
        x=x_col,
        bins=bins,
        kde=False,
        ax=ax,
        edgecolor="black",
        alpha=0.7,
        multiple="stack",
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Frequency")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    save_and_close_plot(output_dir, filename)


def plot_histogram_by_hue(
    data, x_col, hue_col, title, xlabel, output_dir, filename, palette, bins=50, highlight_sample="", sample_values=None
):
    """
    Create and save a histogram plot with hue coloring.

    Args:
        data: Pandas DataFrame containing the data
        x_col: Column name to plot on x-axis
        hue_col: Column name to use for color grouping
        title: Plot title
        xlabel: X-axis label
        output_dir: Directory path where plot will be saved
        filename: Name of output file
        palette: Dictionary mapping hue values to colors
        bins: Number of histogram bins (default: 50)
        highlight_sample: Optional sample ID to highlight with vertical line(s)
        sample_values: Optional list of values for the highlight sample to draw vertical lines

    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(
        data=data,
        x=x_col,
        hue=hue_col,
        bins=bins,
        kde=False,
        palette=palette,
        alpha=0.6,
        edgecolor="black",
        ax=ax,
        multiple="stack",
    )
    
    # Add vertical lines for highlighted sample
    if highlight_sample and sample_values is not None and len(sample_values) > 0:
        for val in sample_values:
            ax.axvline(val, color="magenta", linestyle="--", linewidth=2, alpha=0.8)
        # Add legend entry for highlight sample
        ax.plot([], [], color="magenta", linestyle="--", linewidth=2, label=" " + highlight_sample)
        ax.legend(loc="best", framealpha=0.9)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Frequency")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    save_and_close_plot(output_dir, filename)


def plot_histograms_by_chr_type(df, output_dir, highlight_sample=""):
    """
    Create histograms of all numeric columns colored by chromosome type.

    Args:
        df: Pandas DataFrame with chromosome statistics
        output_dir: Directory path where plots will be saved
        highlight_sample: Optional sample ID to highlight with vertical line(s)

    Returns:
        None
    """
    # Add chromosome type column
    df = df.copy()
    df["chr_type"] = df["chromosome"].apply(get_chromosome_type)

    # Columns to exclude from plotting
    exclude_cols = ["sample", "chromosome", "chr_type"]

    # Get numeric columns
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    plot_cols = [col for col in numeric_cols if col not in exclude_cols]

    # Get highlight sample values if provided
    sample_data = {}
    if highlight_sample and highlight_sample in df["sample"].values:
        highlight_df = df[df["sample"] == highlight_sample]
        for col in plot_cols:
            sample_data[col] = highlight_df[col].tolist()

    # Handle sample_var_map separately - convert to sample_stdev and plot by chr_type
    if "sample_var_map" in plot_cols:
        plot_cols.remove("sample_var_map")
        df["sample_stdev"] = np.sqrt(df["sample_var_map"])
        sample_stdev_vals = None
        if highlight_sample and highlight_sample in df["sample"].values:
            sample_stdev_vals = np.sqrt(highlight_df["sample_var_map"]).tolist()
        plot_histogram_by_hue(
            data=df,
            x_col="sample_stdev",
            hue_col="chr_type",
            title="Distribution of Sample Standard Deviation by Chromosome Type",
            xlabel="Sample Standard Deviation",
            output_dir=output_dir,
            filename="hist_sample_stdev.png",
            palette={"Autosomal": "#1f77b4", "chrX": "#ff7f0e", "chrY": "#2ca02c"},
            bins=50,
            highlight_sample=highlight_sample,
            sample_values=sample_stdev_vals,
        )

    # Color palette for chromosome types
    chr_colors = {"Autosomal": "#1f77b4", "chrX": "#ff7f0e", "chrY": "#2ca02c"}

    # Plot each remaining column
    for col in plot_cols:
        plot_histogram_by_hue(
            data=df,
            x_col=col,
            hue_col="chr_type",
            title=f"Distribution of {format_column_name(col)} by Chromosome Type",
            xlabel=format_column_name(col),
            output_dir=output_dir,
            filename=f"hist_{col}.png",
            palette=chr_colors,
            bins=50,
            highlight_sample=highlight_sample,
            sample_values=sample_data.get(col),
        )


def plot_aneuploid_histograms(df, output_dir, highlight_sample=""):
    """
    Create histograms for aneuploid samples only.

    Args:
        df: Pandas DataFrame with chromosome statistics
        output_dir: Directory path where plots will be saved
        highlight_sample: Optional sample ID to highlight with vertical line(s)

    Returns:
        None
    """
    # Filter to aneuploid samples
    aneuploid_df = df[df["is_aneuploid"]].copy()

    if len(aneuploid_df) == 0:
        print("Warning: No aneuploid samples found, skipping aneuploid histograms")
        return

    # Add chromosome type column for coloring
    aneuploid_df["chr_type"] = aneuploid_df["chromosome"].apply(get_chromosome_type)
    chr_colors = {"Autosomal": "#1f77b4", "chrX": "#ff7f0e", "chrY": "#2ca02c"}

    print(f"\nPlotting aneuploid-only histograms ({len(aneuploid_df)} rows)...")
    
    # Get highlight sample values if in aneuploid data
    sample_data = {}
    if highlight_sample and highlight_sample in aneuploid_df["sample"].values:
        highlight_df = aneuploid_df[aneuploid_df["sample"] == highlight_sample]
        for col in ["copy_number", "mean_cn_probability", "mean_depth", "median_depth", "std_depth"]:
            if col in highlight_df.columns:
                sample_data[col] = highlight_df[col].tolist()

    # Columns to plot
    plot_cols = [
        "chromosome",
        "copy_number",
        "mean_cn_probability",
        "mean_depth",
        "median_depth",
        "std_depth",
    ]

    # Create chromosome histogram
    fig, ax = plt.subplots(figsize=(12, 6))
    # Sort chromosomes properly
    chr_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chr_order = [c for c in chr_order if c in aneuploid_df["chromosome"].unique()]
    sns.countplot(
        data=aneuploid_df,
        x="chromosome",
        order=chr_order,
        ax=ax,
        edgecolor="black",
        alpha=0.7,
    )
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Count")
    ax.set_title("Aneuploid Samples by Chromosome")
    ax.grid(True, alpha=0.3, axis="y")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    diagnostics_dir = os.path.join(output_dir, "diagnostics")
    os.makedirs(diagnostics_dir, exist_ok=True)
    plt.savefig(os.path.join(diagnostics_dir, "hist_aneuploid_chromosome.png"), dpi=300)
    plt.close()
    print("  Created histogram: diagnostics/hist_aneuploid_chromosome.png")

    # Create histograms for numeric columns
    for col in plot_cols[1:]:  # Skip chromosome (already plotted)
        plot_histogram_by_hue(
            data=aneuploid_df,
            x_col=col,
            hue_col="chr_type",
            title=f"Distribution of {format_column_name(col)} by Chromosome Type (Aneuploid Only)",
            xlabel=format_column_name(col),
            output_dir=output_dir,
            filename=f"hist_aneuploid_{col}.png",
            palette=chr_colors,
            bins=30,
            highlight_sample=highlight_sample,
            sample_values=sample_data.get(col),
        )

    # Handle sample_stdev for aneuploid samples
    if "sample_var_map" in aneuploid_df.columns:
        aneuploid_df["sample_stdev"] = np.sqrt(aneuploid_df["sample_var_map"])
        sample_stdev_vals = None
        if highlight_sample and highlight_sample in aneuploid_df["sample"].values:
            highlight_df = aneuploid_df[aneuploid_df["sample"] == highlight_sample]
            sample_stdev_vals = np.sqrt(highlight_df["sample_var_map"]).tolist()
        plot_histogram_by_hue(
            data=aneuploid_df,
            x_col="sample_stdev",
            hue_col="chr_type",
            title="Distribution of Sample Standard Deviation by Chromosome Type (Aneuploid Only)",
            xlabel="Sample Standard Deviation",
            output_dir=output_dir,
            filename="hist_aneuploid_sample_stdev.png",
            palette=chr_colors,
            bins=30,
            highlight_sample=highlight_sample,
            sample_values=sample_stdev_vals,
        )


def load_truth_json(json_path):
    """
    Load truth JSON file with sample to aneuploidy_type mappings.

    Args:
        json_path: Path to JSON file with sample IDs as keys and aneuploidy_types as values

    Returns:
        Dictionary mapping sample IDs to aneuploidy_type names
    """
    with open(json_path, "r") as f:
        return json.load(f)


def assign_sex_and_aneuploidy_types(df, truth_dict):
    """
    Assign sex and aneuploidy types to each sample.

    Args:
        df: Combined DataFrame with all samples and chromosomes
        truth_dict: Dictionary mapping sample IDs to true aneuploidy types (empty dict if not provided)

    Returns:
        DataFrame with columns: sample, sex, chrX_CN, chrY_CN, true_aneuploidy_type, predicted_aneuploidy_type, score
    """
    sample_aneuploidy_types = []

    # Group by sample for efficient iteration
    for sample_id, sample_df in df.groupby("sample"):
        # Get true aneuploidy_type from truth.json
        true_aneuploidy_type = truth_dict.get(sample_id, "NORMAL")

        # Create lookup dictionaries for faster access
        chrom_to_cn = dict(zip(sample_df["chromosome"], sample_df["copy_number"]))
        chrom_to_depth = dict(zip(sample_df["chromosome"], sample_df["median_depth"]))
        chrom_to_aneuploid = dict(
            zip(sample_df["chromosome"], sample_df["is_aneuploid"])
        )

        chrX_cn = chrom_to_cn.get("chrX", 2)
        chrY_cn = chrom_to_cn.get("chrY", 0)
        chrX_depth = chrom_to_depth.get("chrX", 2)
        chrY_depth = chrom_to_depth.get("chrY", 0)
        
        # Determine sex based on chrX and chrY copy numbers
        if chrX_cn == 1 and chrY_cn == 1:
            sex = "MALE"
        elif chrX_cn == 2 and chrY_cn == 0:
            sex = "FEMALE"
        elif chrX_cn == 1 and chrY_cn == 0:
            sex = "TURNER"
        elif chrX_cn == 3 and chrY_cn == 0:
            sex = "TRIPLE_X"
        elif chrX_cn == 2 and chrY_cn == 1:
            sex = "KLINEFELTER"
        elif chrX_cn == 1 and chrY_cn == 2:
            sex = "JACOBS"
        else:
            sex = "OTHER"

        # Predict aneuploidy_type from aneuploidy calls
        aneuploid_chroms = [
            chrom for chrom, is_aneuploid in chrom_to_aneuploid.items() if is_aneuploid
        ]

        if len(aneuploid_chroms) == 0:
            # No aneuploidies detected - assume normal/euploid
            predicted_aneuploidy_type = "NORMAL"
        else:
            # Check if all aneuploidies are on sex chromosomes
            sex_chrom_aneuploidies = [
                c for c in aneuploid_chroms if c in ["chrX", "chrY"]
            ]
            autosomal_aneuploidies = [
                c for c in aneuploid_chroms if c not in ["chrX", "chrY"]
            ]

            # If there are multiple autosomal aneuploidies, or mix of autosomal and sex, call it MULTIPLE
            if len(autosomal_aneuploidies) > 1 or (
                len(autosomal_aneuploidies) > 0 and len(sex_chrom_aneuploidies) > 0
            ):
                predicted_aneuploidy_type = "MULTIPLE"
            # If only sex chromosome aneuploidies, determine aneuploidy_type from chrX/chrY counts
            elif len(sex_chrom_aneuploidies) > 0 and len(autosomal_aneuploidies) == 0:
                if chrX_cn == 2 and chrY_cn == 1:
                    predicted_aneuploidy_type = "KLINEFELTER"  # XXY
                elif chrX_cn == 3 and chrY_cn == 0:
                    predicted_aneuploidy_type = "TRIPLE_X"  # XXX
                elif chrX_cn == 1 and chrY_cn == 0:
                    predicted_aneuploidy_type = "TURNER"  # X0
                elif chrX_cn == 1 and chrY_cn == 2:
                    predicted_aneuploidy_type = "JACOBS"  # XYY
                else:
                    predicted_aneuploidy_type = "OTHER"
            # Single autosomal aneuploidy
            elif len(autosomal_aneuploidies) == 1:
                chrom = autosomal_aneuploidies[0]
                cn = chrom_to_cn[chrom]

                if chrom == "chr13" and cn == 3:
                    predicted_aneuploidy_type = "TRISOMY_13"
                elif chrom == "chr18" and cn == 3:
                    predicted_aneuploidy_type = "TRISOMY_18"
                elif chrom == "chr21" and cn == 3:
                    predicted_aneuploidy_type = "TRISOMY_21"
                else:
                    predicted_aneuploidy_type = "OTHER"
            else:
                predicted_aneuploidy_type = "OTHER"

        # Calculate minimum mean_cn_probability across all chromosomes
        min_cn_prob = sample_df["mean_cn_probability"].min()

        sample_aneuploidy_types.append(
            {
                "sample": sample_id,
                "sex": sex,
                "chrX_CN": chrX_cn,
                "chrY_CN": chrY_cn,
                "chrX_depth": chrX_depth,
                "chrY_depth": chrY_depth,
                "true_aneuploidy_type": true_aneuploidy_type,
                "predicted_aneuploidy_type": predicted_aneuploidy_type,
                "score": min_cn_prob,
            }
        )

    return pd.DataFrame(sample_aneuploidy_types)


def calculate_metrics(aneuploidy_type_df, output_dir):
    """
    Calculate and report sensitivity, specificity, and precision metrics.

    Args:
        aneuploidy_type_df: DataFrame with true_aneuploidy_type and predicted_aneuploidy_type columns
        output_dir: Directory to save metrics report

    Returns:
        None
    """
    print("\n" + "=" * 80)
    print("aneuploidy_type PREDICTION METRICS")
    print("=" * 80)

    # Overall accuracy
    correct = (
        aneuploidy_type_df["true_aneuploidy_type"] == aneuploidy_type_df["predicted_aneuploidy_type"]
    ).sum()
    total = len(aneuploidy_type_df)
    accuracy = correct / total if total > 0 else 0
    print(
        f"\nOverall Accuracy: {correct}/{total} = {accuracy:.4f} ({100 * accuracy:.2f}%)"
    )

    # Per-aneuploidy_type metrics
    print("\n" + "-" * 80)
    print("Per-aneuploidy_type Classification Report:")
    print("-" * 80)

    # Get classification report
    classification_report(
        aneuploidy_type_df["true_aneuploidy_type"],
        aneuploidy_type_df["predicted_aneuploidy_type"],
        zero_division=0,
        output_dict=True,
    )

    # Print formatted report
    report_str = classification_report(
        aneuploidy_type_df["true_aneuploidy_type"],
        aneuploidy_type_df["predicted_aneuploidy_type"],
        zero_division=0,
    )
    print(report_str)

    # Create confusion matrix
    labels = sorted(aneuploidy_type_df["true_aneuploidy_type"].unique())
    cm = confusion_matrix(
        aneuploidy_type_df["true_aneuploidy_type"],
        aneuploidy_type_df["predicted_aneuploidy_type"],
        labels=labels,
    )

    print("\n" + "-" * 80)
    print("Confusion Matrix:")
    print("-" * 80)
    cm_df = pd.DataFrame(cm, index=labels, columns=labels)
    print(cm_df)

    # Save detailed results
    metrics_path = os.path.join(output_dir, "aneuploidy_type_predictions.tsv")
    aneuploidy_type_df.to_csv(metrics_path, sep="\t", index=False)
    print(f"\nDetailed predictions saved to: {metrics_path}")

    # Save metrics report
    report_path = os.path.join(output_dir, "metrics_report.txt")
    with open(report_path, "w") as f:
        f.write("=" * 80 + "\n")
        f.write("aneuploidy_type PREDICTION METRICS\n")
        f.write("=" * 80 + "\n\n")
        f.write(
            f"Overall Accuracy: {correct}/{total} = {accuracy:.4f} ({100 * accuracy:.2f}%)\n\n"
        )
        f.write("-" * 80 + "\n")
        f.write("Per-aneuploidy_type Classification Report:\n")
        f.write("-" * 80 + "\n")
        f.write(report_str)
        f.write("\n" + "-" * 80 + "\n")
        f.write("Confusion Matrix:\n")
        f.write("-" * 80 + "\n")
        f.write(cm_df.to_string())

    print(f"Metrics report saved to: {report_path}")
    print("=" * 80 + "\n")


def load_exclusion_ids(file_path):
    """
    Load excluded case IDs from text file.

    Args:
        file_path: Path to text file with exclusion IDs (one per line)

    Returns:
        List of exclusion ID strings
    """
    with open(file_path, "r") as f:
        exclusion_ids = [line.strip() for line in f if line.strip()]
        return exclusion_ids


def export_aneuploid_data(df, output_path):
    """
    Export rows where is_aneuploid=True to a TSV file.

    Args:
        df: DataFrame with aneuploidy data
        output_path: Path where TSV file will be saved

    Returns:
        None

    Note:
        Converts sample_var_map to inferred_sample_std and removes sample_var_map and chr_type columns
    """
    aneuploid_df = df[df["is_aneuploid"]].copy()
    aneuploid_df["inferred_sample_std"] = np.sqrt(aneuploid_df["sample_var_map"])
    aneuploid_df = aneuploid_df.drop(columns=["sample_var_map"], errors="ignore")
    aneuploid_df = aneuploid_df.drop(columns=["chr_type"], errors="ignore")
    aneuploid_df.to_csv(output_path, sep="\t", index=False)
    print(f"\nExported {len(aneuploid_df)} aneuploid rows to {output_path}")


def generate_median_depth_distributions(df, output_dir, highlight_sample=""):
    """
    Generate median depth distribution statistics and plots with MAD error bars.

    Args:
        df: DataFrame with chromosome statistics including median_depth and mad_depth columns
        output_dir: Directory to save outputs
        highlight_sample: Optional sample ID to highlight in plots

    Returns:
        None
    """
    print("\nGenerating median depth distribution plots...")

    # Check if required columns exist
    if "median_depth" not in df.columns or "mad_depth" not in df.columns:
        print("Warning: 'median_depth' or 'mad_depth' column not found, skipping CN denoising outputs")
        return

    # Create stats dataframe with required columns
    stats_df = df[["sample", "chromosome", "copy_number", "mean_depth", "median_depth", "mad_depth"]].copy()
    
    # Save stats to TSV
    stats_path = os.path.join(output_dir, "median_depth_distribution_stats.tsv")
    stats_df.to_csv(stats_path, sep="\t", index=False)
    print(f"  Saved median depth distribution stats to {stats_path}")

    # Generate error bar plots (one per chromosome)
    unique_chromosomes = sorted(stats_df["chromosome"].unique(), 
                                key=lambda x: (x in ["chrX", "chrY"], 
                                             int(x.replace("chr", "")) if x.replace("chr", "").isdigit() else 0,
                                             x))
    
    pdf_path = os.path.join(output_dir, "median_depth_distributions.pdf")
    with PdfPages(pdf_path) as pdf:
        for chrom in unique_chromosomes:
            chrom_data = stats_df[stats_df["chromosome"] == chrom].copy()
            chrom_data = chrom_data.sort_values("median_depth")

            fig, ax = plt.subplots(figsize=(14, 6))
            
            # Find highlight sample index if present
            highlight_idx = None
            if highlight_sample and highlight_sample in chrom_data["sample"].values:
                highlight_idx = chrom_data[chrom_data["sample"] == highlight_sample].index[0]
                highlight_pos = chrom_data.index.get_loc(highlight_idx)
            
            # Plot all samples
            ax.errorbar(
                x=range(len(chrom_data)),
                y=chrom_data["median_depth"],
                yerr=chrom_data["mad_depth"],
                fmt="o",
                markerfacecolor="blue",
                markeredgecolor="blue",
                ecolor="black",
                capsize=2,
                linestyle="none",
                alpha=0.7
            )
            
            # Highlight the specified sample if present
            if highlight_idx is not None:
                highlight_data = chrom_data.loc[highlight_idx]
                ax.errorbar(
                    x=[highlight_pos],
                    y=[highlight_data["median_depth"]],
                    yerr=[highlight_data["mad_depth"]],
                    fmt="^",
                    markerfacecolor="magenta",
                    markeredgecolor="black",
                    ecolor="magenta",
                    capsize=4,
                    markersize=10,
                    linewidth=2,
                    linestyle="none",
                    label=" " + highlight_sample
                )
                ax.legend(loc="best", framealpha=0.9)

            ax.set_xlabel("Sample Index")
            ax.set_ylabel("Median Depth")
            ax.set_title(f"Median Depth for {chrom} with MAD as Error Bars")
            ax.set_xticks(range(len(chrom_data)))
            ax.set_xticklabels(chrom_data["sample"], rotation=90, fontsize=6)
            plt.tight_layout()
            pdf.savefig()
            plt.close()
    
    print(f"  Saved median depth distribution plots to {pdf_path}")


def plot_training_loss(training_loss_df, output_dir):
    """
    Plot training loss over epochs.
    
    Args:
        training_loss_df: DataFrame with 'epoch' and 'elbo' columns
        output_dir: Directory to save plot
        
    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(training_loss_df["epoch"], training_loss_df["elbo"])
    ax.set_xlabel("Epoch")
    ax.set_ylabel("ELBO")
    ax.set_title("Training Loss")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    diagnostics_dir = os.path.join(output_dir, "diagnostics")
    os.makedirs(diagnostics_dir, exist_ok=True)
    output_path = os.path.join(diagnostics_dir, "training_loss.png")
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Created plot: diagnostics/training_loss.png")


def add_chromosome_labels(ax, chr_array, x_transformed=None):
    """
    Add chromosome labels to x-axis of a plot.
    
    Args:
        ax: matplotlib axis
        chr_array: array of chromosome names for each bin
        x_transformed: optional array of transformed x coordinates (for equal-width chromosomes)
        
    Returns:
        None
    """
    # Find chromosome boundaries
    chr_changes = np.where(chr_array[:-1] != chr_array[1:])[0] + 1
    chr_boundaries = np.concatenate([[0], chr_changes, [len(chr_array)]])
    
    # Calculate midpoint of each chromosome for labels
    chr_labels = []
    chr_positions = []
    chr_boundary_positions = []
    
    for i in range(len(chr_boundaries) - 1):
        start = chr_boundaries[i]
        end = chr_boundaries[i + 1]
        
        if x_transformed is not None:
            # Use transformed coordinates - chromosome centers are at i+0.5
            chr_positions.append(i + 0.5)
            if i > 0:
                chr_boundary_positions.append(i)
        else:
            # Use original bin indices
            mid = (start + end) / 2
            chr_positions.append(mid)
            if start > 0:
                chr_boundary_positions.append(start)
        
        chr_name = chr_array[start]
        chr_labels.append(str(chr_name).replace("chr", ""))
    
    # Add vertical lines at chromosome boundaries
    for boundary in chr_boundary_positions:
        ax.axvline(boundary, color="gray", linestyle="--", alpha=1, linewidth=1)
    
    # Set x-axis labels
    ax.set_xticks(chr_positions)
    ax.set_xticklabels(chr_labels, rotation=0, ha="center")
    ax.set_xlabel("Chromosome")


def plot_sample(sample_data, output_dir, aneuploid_chrs=None):
    """
    Combined plot showing observed depth, CN posterior, and sample variance.
    
    Args:
        sample_data: DataFrame with bin data for a single sample
        output_dir: Directory to save plots
        aneuploid_chrs: List of tuples (chr_name, cn, prob) for aneuploid chromosomes
        
    Returns:
        None
    """
    sample_name = sample_data["sample"].iloc[0]
    is_aneuploid = aneuploid_chrs is not None and len(aneuploid_chrs) > 0
    
    fig, axes = plt.subplots(3, 1, figsize=(8, 8))
    
    # Add sample name as figure title
    fig.text(
        0.1 / 8,
        0.98,
        f"{sample_name}",
        fontsize=12,
        fontweight="bold",
        va="top",
        ha="left",
    )
    
    # Extract data
    observed = sample_data["observed_depth"].values
    cn_map = sample_data["cn_map"].values
    chr_array = sample_data["chr"].values
    
    # Get CN probabilities
    cn_probs = sample_data[["cn_prob_0", "cn_prob_1", "cn_prob_2", "cn_prob_3", "cn_prob_4", "cn_prob_5"]].values
    
    # Get sample variance (same for all bins of this sample)
    sample_var = sample_data["sample_var"].iloc[0]
    
    # Colors for CN states
    colors = ["#004D40", "#FFC107", "#1E88E5", "#D81B60", "#38006B", "#FF6D00"]
    
    # Create normalized x-axis where each chromosome has equal width
    chr_changes = np.where(chr_array[:-1] != chr_array[1:])[0] + 1
    chr_boundaries = np.concatenate([[0], chr_changes, [len(chr_array)]])
    
    x_transformed = np.zeros(len(observed))
    for i in range(len(chr_boundaries) - 1):
        start_idx = chr_boundaries[i]
        end_idx = chr_boundaries[i + 1]
        n_bins_in_chr = end_idx - start_idx
        x_transformed[start_idx:end_idx] = i + np.linspace(0, 1, n_bins_in_chr, endpoint=False)
    
    # Plot 1: Observed depth with MAP copy number calls overlaid
    ax = axes[0]
    ax.plot(x_transformed, cn_map, "o", alpha=0.5, markersize=4, 
            label="Predicted bin copy number", color="black")
    ax.plot(x_transformed, observed, "r-", linewidth=1, alpha=0.7, 
            label="Normalized read depth")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=2, borderaxespad=0)
    ax.grid(True, axis="y", alpha=1, linestyle="-", linewidth=1)
    ax.set_ylim([-0.5, 5.5])
    ax.set_xlim([x_transformed.min(), x_transformed.max()])
    
    # Plot 2: Stacked area plot of CN probabilities
    ax = axes[1]
    ax.stackplot(x_transformed, cn_probs.T, 
                 labels=[f"CN={i}" for i in range(6)], alpha=0.7, colors=colors)
    ax.set_ylabel("Copy Number Probability")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=6, borderaxespad=0)
    ax.set_ylim([0, 1])
    ax.set_xlim([x_transformed.min(), x_transformed.max()])
    
    # Plot 3: Sample variance distribution (placeholder - we'll need all samples' variance)
    # For now, just show this sample's variance
    ax = axes[2]
    ax.axvline(np.sqrt(sample_var), color="red", linestyle="--", linewidth=2,
               label=f"This sample: {np.sqrt(sample_var):.3f}")
    ax.set_xlabel("Sample stdev")
    ax.set_ylabel("Count")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.02), borderaxespad=0)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, np.sqrt(sample_var) * 2])
    
    # Highlight aneuploid chromosomes in first 2 panels
    if is_aneuploid and aneuploid_chrs:
        aneuploid_chr_names = set([chr_name for chr_name, _, _ in aneuploid_chrs])
        
        for chr_name in aneuploid_chr_names:
            chr_mask = chr_array == chr_name
            chr_indices = np.where(chr_mask)[0]
            
            if len(chr_indices) > 0:
                chr_start_transformed = x_transformed[chr_indices[0]]
                chr_end_transformed = x_transformed[chr_indices[-1]]
                
                for ax in axes[:2]:
                    ax.axvspan(chr_start_transformed, chr_end_transformed, 
                              alpha=0.15, color="red", zorder=0)
    
    # Add chromosome labels to first 2 panels
    for ax in axes[:2]:
        add_chromosome_labels(ax, chr_array, x_transformed=x_transformed)
    
    plt.tight_layout()
    sample_name_clean = sample_name.replace("/", "_").replace(" ", "_")
    filename = f"{sample_name_clean}.png"
    
    subdir = "sample_plots"
    subdir_path = os.path.join(output_dir, subdir)
    os.makedirs(subdir_path, exist_ok=True)
    
    output_path = os.path.join(subdir_path, filename)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_sample_with_variance_distribution(sample_data, all_sample_vars, output_dir, aneuploid_chrs=None):
    """
    Combined plot showing observed depth, CN posterior, and sample variance distribution.
    
    Args:
        sample_data: DataFrame with bin data for a single sample
        all_sample_vars: Array of all samples' variance values
        output_dir: Directory to save plots
        aneuploid_chrs: List of tuples (chr_name, cn, prob) for aneuploid chromosomes
        
    Returns:
        None
    """
    sample_name = sample_data["sample"].iloc[0]
    is_aneuploid = aneuploid_chrs is not None and len(aneuploid_chrs) > 0
    
    fig, axes = plt.subplots(3, 1, figsize=(8, 8))
    
    # Add sample name as figure title
    fig.text(
        0.1 / 8,
        0.98,
        f"{sample_name}",
        fontsize=12,
        fontweight="bold",
        va="top",
        ha="left",
    )
    
    # Extract data
    observed = sample_data["observed_depth"].values
    cn_map = sample_data["cn_map"].values
    chr_array = sample_data["chr"].values
    
    # Get CN probabilities
    cn_probs = sample_data[["cn_prob_0", "cn_prob_1", "cn_prob_2", "cn_prob_3", "cn_prob_4", "cn_prob_5"]].values
    
    # Get sample variance
    sample_var = sample_data["sample_var"].iloc[0]
    
    # Colors for CN states
    colors = ["#004D40", "#FFC107", "#1E88E5", "#D81B60", "#38006B", "#FF6D00"]
    
    # Create normalized x-axis where each chromosome has equal width
    chr_changes = np.where(chr_array[:-1] != chr_array[1:])[0] + 1
    chr_boundaries = np.concatenate([[0], chr_changes, [len(chr_array)]])
    
    x_transformed = np.zeros(len(observed))
    for i in range(len(chr_boundaries) - 1):
        start_idx = chr_boundaries[i]
        end_idx = chr_boundaries[i + 1]
        n_bins_in_chr = end_idx - start_idx
        x_transformed[start_idx:end_idx] = i + np.linspace(0, 1, n_bins_in_chr, endpoint=False)
    
    # Plot 1: Observed depth with MAP copy number calls overlaid
    ax = axes[0]
    ax.plot(x_transformed, cn_map, "o", alpha=0.5, markersize=4, 
            label="Predicted bin copy number", color="black")
    ax.plot(x_transformed, observed, "r-", linewidth=1, alpha=0.7, 
            label="Normalized read depth")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=2, borderaxespad=0)
    ax.grid(True, axis="y", alpha=1, linestyle="-", linewidth=1)
    ax.set_ylim([-0.5, 5.5])
    ax.set_xlim([x_transformed.min(), x_transformed.max()])
    
    # Plot 2: Stacked area plot of CN probabilities
    ax = axes[1]
    ax.stackplot(x_transformed, cn_probs.T, 
                 labels=[f"CN={i}" for i in range(6)], alpha=0.7, colors=colors)
    ax.set_ylabel("Copy Number Probability")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=6, borderaxespad=0)
    ax.set_ylim([0, 1])
    ax.set_xlim([x_transformed.min(), x_transformed.max()])
    
    # Plot 3: Sample variance distribution
    ax = axes[2]
    ax.hist(np.sqrt(all_sample_vars), bins=30, alpha=0.7, edgecolor="black", 
            linewidth=0.5, color="gray")
    ax.axvline(np.sqrt(sample_var), color="red", linestyle="--", linewidth=2,
               label=f"This sample: {np.sqrt(sample_var):.3f}")
    ax.set_xlabel("Sample stdev")
    ax.set_ylabel("Count")
    ax.set_yscale("log")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.02), borderaxespad=0)
    ax.grid(True, alpha=0.3)
    
    # Highlight aneuploid chromosomes in first 2 panels
    if is_aneuploid and aneuploid_chrs:
        aneuploid_chr_names = set([chr_name for chr_name, _, _ in aneuploid_chrs])
        
        for chr_name in aneuploid_chr_names:
            chr_mask = chr_array == chr_name
            chr_indices = np.where(chr_mask)[0]
            
            if len(chr_indices) > 0:
                chr_start_transformed = x_transformed[chr_indices[0]]
                chr_end_transformed = x_transformed[chr_indices[-1]]
                
                for ax in axes[:2]:
                    ax.axvspan(chr_start_transformed, chr_end_transformed, 
                              alpha=0.15, color="red", zorder=0)
    
    # Add chromosome labels to first 2 panels
    for ax in axes[:2]:
        add_chromosome_labels(ax, chr_array, x_transformed=x_transformed)
    
    plt.tight_layout()
    sample_name_clean = sample_name.replace("/", "_").replace(" ", "_")
    aneu_suffix = "ANEU" if is_aneuploid else "NORMAL"
    filename = f"sample_{aneu_suffix}_{sample_name_clean}.png"
    
    subdir = "sample_plots"
    subdir_path = os.path.join(output_dir, subdir)
    os.makedirs(subdir_path, exist_ok=True)
    
    output_path = os.path.join(subdir_path, filename)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_bin_variance_bias(bin_data, output_dir):
    """
    Plot bin bias and bin variance posteriors (global across all samples).
    
    Args:
        bin_data: DataFrame with bin-level data including bin_bias and bin_var columns
        output_dir: Directory to save plot
        
    Returns:
        None
    """
    # Get unique bins (one row per bin, across all samples)
    # Group by bin and take first occurrence
    unique_bins = bin_data.groupby(["chr", "start", "end"]).first().reset_index()
    unique_bins = unique_bins.sort_values(["chr", "start"])
    
    fig, axes = plt.subplots(2, 1, figsize=(16, 8))
    
    bin_bias = unique_bins["bin_bias"].values
    bin_var = unique_bins["bin_var"].values
    chr_array = unique_bins["chr"].values
    
    x = np.arange(len(bin_bias))
    
    # Plot 1: Bin bias
    ax = axes[0]
    ax.plot(x, bin_bias, "o", alpha=0.5, markersize=3, color="black", label="Bin bias")
    ax.axhline(1.0, color="red", linestyle="--", alpha=0.5, label="No bias")
    ax.set_ylabel("Bin mean bias")
    ax.set_title("Bin Posteriors")
    ax.set_xlim([x.min(), x.max()])
    
    # Plot 2: Bin variance
    ax = axes[1]
    ax.plot(x, np.sqrt(bin_var), "o", alpha=0.5, markersize=3, color="purple", label="Bin stdev")
    ax.set_ylabel("Bin stdev")
    ax.set_xlim([x.min(), x.max()])
    
    # Add chromosome labels to both panels
    for ax in axes:
        add_chromosome_labels(ax, chr_array)
    
    plt.tight_layout()
    diagnostics_dir = os.path.join(output_dir, "diagnostics")
    os.makedirs(diagnostics_dir, exist_ok=True)
    filename = "bin_posteriors.png"
    output_path = os.path.join(diagnostics_dir, filename)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Created plot: diagnostics/bin_posteriors.png")


def plot_sex_assignments(sex_df, output_dir, highlight_sample=""):
    """
    Create scatter plot of chrX vs chrY copy number with sex assignments.
    Replicates estimatePloidy.R sex_assignments.png plot.
    
    Args:
        sex_df: DataFrame with columns: sample, sex, chrX_CN, chrY_CN
        output_dir: Directory to save plots
        highlight_sample: Optional sample ID to highlight
        
    Returns:
        None
    """
    print("\nGenerating sex assignment scatter plot...")
    
    # Use sex column directly for assignments
    sex_df = sex_df.copy()
    
    # Define colors matching R script
    sex_colors = {
        "MALE": "#00BFF4",
        "FEMALE": "#fd8eff",
        "TURNER": "#e02006",
        "TRIPLE_X": "#7B2AB3",
        "KLINEFELTER": "#FF6A09",
        "JACOBS": "#29840f",
        "OTHER": "#8F1336"
    }
    
    # Create plot
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Set axis limits
    ax_lim = 3
    ax.set_xlim(-0.1, ax_lim + 0.1)
    ax.set_ylim(-0.1, ax_lim + 0.1)
    
    # Add grid lines
    for i in range(ax_lim + 1):
        ax.axhline(i, color="gray", linestyle=":", alpha=0.5)
        ax.axvline(i, color="gray", linestyle=":", alpha=0.5)
    
    # Add karyotype labels in background
    for x_cn in range(ax_lim + 1):
        for y_cn in range(ax_lim + 1):
            karyo = "X" * x_cn + "Y" * y_cn
            ax.text(x_cn, y_cn, karyo, fontsize=16, fontweight="bold",
                   color="lightgray", ha="center", va="center", alpha=0.8)
    
    # Plot points by sex assignment
    for assignment in sex_df["sex"].unique():
        subset = sex_df[sex_df["sex"] == assignment]
        is_highlight = subset["sample"] == highlight_sample
        
        # Plot non-highlighted samples
        if not is_highlight.any() or len(subset) > 1:
            non_highlight = subset[subset["sample"] != highlight_sample]
            if len(non_highlight) > 0:
                ax.scatter(non_highlight["chrX_depth"], non_highlight["chrY_depth"],
                          c=sex_colors.get(assignment, "#8F1336"),
                          s=50, alpha=0.7, label=assignment, edgecolors="none")
        
        # Plot highlighted sample separately
        if is_highlight.any():
            highlight_row = subset[is_highlight]
            ax.scatter(highlight_row["chrX_depth"], highlight_row["chrY_depth"],
                      c=sex_colors.get(assignment, "#8F1336"),
                      s=100, alpha=1.0, marker="^", edgecolors="black", linewidths=2)
    
    # Add crosshairs for highlighted sample
    if highlight_sample and highlight_sample in sex_df["sample"].values:
        highlight_row = sex_df[sex_df["sample"] == highlight_sample].iloc[0]
        ax.axhline(highlight_row["chrY_depth"], color="gray", linestyle="-", linewidth=2, alpha=0.7, label=" " + highlight_sample)
        ax.axvline(highlight_row["chrX_depth"], color="gray", linestyle="-", linewidth=2, alpha=0.7)
        # Show both sex assignment legend and highlight sample
        ax.legend(loc="upper right", framealpha=0.9)
    else:
        # Only show sex assignment legend
        ax.legend(loc="upper right", framealpha=0.9)
    
    ax.set_xlabel("chrX Normalized Depth", fontsize=12)
    ax.set_ylabel("chrY Normalized Depth", fontsize=12)
    
    plt.tight_layout()
    output_path = os.path.join(output_dir, "sex_assignments.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Created plot: sex_assignments.png")


def plot_cn_per_contig_boxplot(df, output_dir, sex_subset=None, connect_samples=True, highlight_sample=""):
    """
    Create boxplot with jitter showing copy number per contig.
    Replicates estimatePloidy.R estimated_CN_per_contig plots.
    
    Args:
        df: DataFrame with chromosome statistics
        output_dir: Directory to save plots
        sex_subset: "male", "female", "all", or None (for all samples)
        connect_samples: Whether to connect points across chromosomes with lines
        highlight_sample: Optional sample ID to highlight
        
    Returns:
        None
    """
    # Prepare data
    plot_df = df.copy()
    
    # Filter by sex if specified
    if sex_subset == "male":
        # Males: chrX_CN < 2
        male_samples = []
        for sample, sample_df in df.groupby("sample"):
            chrX_cn = sample_df[sample_df["chromosome"] == "chrX"]["copy_number"].values
            if len(chrX_cn) > 0 and chrX_cn[0] < 2:
                male_samples.append(sample)
        plot_df = plot_df[plot_df["sample"].isin(male_samples)]
        suffix = "chrX_lessThan_2copies"
    elif sex_subset == "female":
        # Females: chrX_CN >= 2
        female_samples = []
        for sample, sample_df in df.groupby("sample"):
            chrX_cn = sample_df[sample_df["chromosome"] == "chrX"]["copy_number"].values
            if len(chrX_cn) > 0 and chrX_cn[0] >= 2:
                female_samples.append(sample)
        plot_df = plot_df[plot_df["sample"].isin(female_samples)]
        suffix = "chrX_atLeast_2copies"
    else:
        suffix = "all_samples"
    
    # Create chromosome order
    chr_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chr_order = [c for c in chr_order if c in plot_df["chromosome"].unique()]
    
    # Pivot to wide format for plotting
    plot_wide = plot_df.pivot(index="sample", columns="chromosome", values="median_depth")
    plot_wide = plot_wide[chr_order]
    
    # Also pivot copy_number for color determination
    cn_wide = plot_df.pivot(index="sample", columns="chromosome", values="copy_number")
    cn_wide = cn_wide[chr_order]
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Set y-axis limits
    y_max = max(4, plot_df["median_depth"].max())
    ax.set_ylim(0, y_max)
    
    # Add shading for alternating chromosomes
    for i in range(0, len(chr_order), 2):
        ax.axvspan(i - 0.5, i + 0.5, color="lightgray", alpha=0.3, zorder=0)
    
    # Add horizontal gridlines
    for y in np.arange(0, y_max + 0.5, 0.5):
        if y == int(y):
            ax.axhline(y, color="gray", linestyle="-", alpha=0.5, linewidth=0.8, zorder=1)
        else:
            ax.axhline(y, color="gray", linestyle=":", alpha=0.3, linewidth=0.5, zorder=1)
    
    # Add reference line at ploidy=2
    ax.axhline(2, color="black", linestyle="-", linewidth=1.5, zorder=2)
    
    # Connect samples with lines if requested
    if connect_samples:
        for sample in plot_wide.index:
            sample_vals = plot_wide.loc[sample].values
            x_coords = np.arange(len(chr_order))
            
            # Only plot line for non-highlighted samples or use gray for highlighted
            if sample == highlight_sample:
                continue
            else:
                ax.plot(x_coords, sample_vals, color="gray", alpha=0.2, linewidth=0.4, zorder=3)
    
    # Plot jittered points with color based on aneuploidy
    for i, chrom in enumerate(chr_order):
        chrom_plot_df = plot_df[plot_df["chromosome"] == chrom].copy()
        chrom_data = chrom_plot_df["median_depth"].values
        chrom_cn = chrom_plot_df["copy_number"].values
        
        # Add jitter
        x_jitter = np.random.normal(i, 0.1, size=len(chrom_data))
        
        # Determine expected ploidy range for this chromosome
        if chrom.startswith("chr") and chrom[3:].isdigit():
            # Autosome: expected = 2
            expected_min = 2
            expected_max = 2
        elif chrom == "chrX":
            # chrX: expected range [1, 2]
            expected_min = 1
            expected_max = 2
        elif chrom == "chrY":
            # chrY: expected range [0, 1]
            expected_min = 0
            expected_max = 1
        else:
            # Default to diploid
            expected_min = 2
            expected_max = 2
        
        # Assign colors based on copy number
        colors = []
        for cn in chrom_cn:
            if cn > expected_max:
                colors.append("blue")  # Gain
            elif cn < expected_min:
                colors.append("red")   # Loss
            else:
                colors.append("#838393")  # Normal
        
        ax.scatter(x_jitter, chrom_data, s=10, alpha=0.5, c=colors, zorder=4)
    
    # Overlay boxplots
    boxplot_data = [plot_df[plot_df["chromosome"] == chrom]["median_depth"].values for chrom in chr_order]
    bp = ax.boxplot(boxplot_data, positions=np.arange(len(chr_order)), widths=0.6,
                    patch_artist=False, showcaps=True, showfliers=False,
                    boxprops=dict(linewidth=1.5, color="black"),
                    whiskerprops=dict(linewidth=1.5, color="black"),
                    medianprops=dict(linewidth=1.5, color="black"),
                    capprops=dict(linewidth=0),
                    zorder=5)
    
    # Highlight specific sample if requested
    if highlight_sample and highlight_sample in plot_wide.index:
        sample_vals = plot_wide.loc[highlight_sample].values
        ax.plot(np.arange(len(chr_order)), sample_vals, color="magenta", linewidth=2, zorder=6, label=" " + highlight_sample)
        ax.legend(loc="upper right", framealpha=0.9)
    
    # Format axes
    ax.set_xticks(np.arange(len(chr_order)))
    ax.set_xticklabels([c.replace("chr", "") for c in chr_order], rotation=90)
    ax.set_xlabel("Chromosome", fontsize=12)
    ax.set_ylabel("Normalized Depth", fontsize=12)
    
    # Add sample count text
    n_samples = plot_wide.shape[0]
    n_incomplete = plot_wide.isnull().any(axis=1).sum()
    ax.text(0.02, 0.98, f"N={n_samples} samples ({n_incomplete} incomplete)",
           transform=ax.transAxes, va="top", ha="left", fontsize=10,
           bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    
    plt.tight_layout()
    
    # Save with appropriate filename
    contour_suffix = "with_contours" if connect_samples else "no_contours"
    filename = f"estimated_CN_per_contig.{suffix}.{contour_suffix}.png"
    raw_depth_dir = os.path.join(output_dir, "raw_depth_contig")
    os.makedirs(raw_depth_dir, exist_ok=True)
    output_path = os.path.join(raw_depth_dir, filename)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Created plot: {raw_depth_dir}/{filename}")


def plot_cn_per_bin_chromosome(bin_stats_df, output_dir, chromosome, sex_subset=None, highlight_sample=""):
    """
    Create boxplot with jitter showing copy number per bin for a specific chromosome.
    Replicates estimatePloidy.R estimated_CN_per_bin plots.
    
    Args:
        bin_stats_df: DataFrame with bin-level statistics
        output_dir: Directory to save plots
        chromosome: Chromosome to plot (e.g., "chr1", "chrX")
        sex_subset: "male", "female", or None
        highlight_sample: Optional sample ID to highlight
        
    Returns:
        None
    """
    # Filter to chromosome
    chr_df = bin_stats_df[bin_stats_df["chr"] == chromosome].copy()
    
    if len(chr_df) == 0:
        return
    
    # Determine which column to plot
    # Sex chromosomes (chrX, chrY) always use observed_depth
    # Autosomes use observed_depth for all_samples, cn_map for male/female subsets
    if chromosome in ["chrX", "chrY"]:
        plot_col = "observed_depth"
    elif sex_subset is None:
        plot_col = "observed_depth"
    else:
        plot_col = "cn_map"
    
    # Filter by sex if specified
    if sex_subset == "male":
        # Get male samples (chrX < 2 copies)
        chr_df_wide = bin_stats_df.pivot_table(index="sample", columns="chr", values="cn_map", aggfunc="mean")
        if "chrX" in chr_df_wide.columns:
            male_samples = chr_df_wide[chr_df_wide["chrX"] < 2].index.tolist()
            chr_df = chr_df[chr_df["sample"].isin(male_samples)]
        suffix = "chrX_lessThan_2copies"
    elif sex_subset == "female":
        # Get female samples (chrX >= 2 copies)
        chr_df_wide = bin_stats_df.pivot_table(index="sample", columns="chr", values="cn_map", aggfunc="mean")
        if "chrX" in chr_df_wide.columns:
            female_samples = chr_df_wide[chr_df_wide["chrX"] >= 2].index.tolist()
            chr_df = chr_df[chr_df["sample"].isin(female_samples)]
        suffix = "chrX_atLeast_2copies"
    else:
        suffix = "all_samples"
    
    # Create bin position index
    chr_df = chr_df.sort_values("start")
    bins = chr_df[["start", "end"]].drop_duplicates().sort_values("start").reset_index(drop=True)
    bins["bin_idx"] = bins.index
    chr_df = chr_df.merge(bins[["start", "end", "bin_idx"]], on=["start", "end"], how="left")
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Set y-axis limits
    y_max = 5
    ax.set_ylim(0, y_max)
    
    # Add horizontal gridlines
    for y in np.arange(0, y_max + 0.5, 0.5):
        if y == int(y):
            ax.axhline(y, color="gray", linestyle="-", alpha=0.5, linewidth=0.8)
        else:
            ax.axhline(y, color="gray", linestyle=":", alpha=0.3, linewidth=0.5)
    
    # Add reference line at ploidy=2
    ax.axhline(2, color="black", linestyle="-", linewidth=1.5)
    
    # Connect samples with lines
    for sample in chr_df["sample"].unique():
        sample_df = chr_df[chr_df["sample"] == sample].sort_values("bin_idx")
        if sample == highlight_sample:
            continue
        ax.plot(sample_df["bin_idx"], sample_df[plot_col], color="gray", alpha=0.2, linewidth=2)
    
    # Plot jittered points for each bin
    bin_indices = sorted(chr_df["bin_idx"].unique())
    for bin_idx in bin_indices:
        bin_data = chr_df[chr_df["bin_idx"] == bin_idx][plot_col].values
        
        # Add jitter
        x_jitter = np.random.normal(bin_idx, 0.1, size=len(bin_data))
        ax.scatter(x_jitter, bin_data, s=10, alpha=0.5, color="#838393")
    
    # Overlay boxplots - only for bins with data
    boxplot_data = []
    boxplot_positions = []
    for idx in bin_indices:
        bin_data = chr_df[chr_df["bin_idx"] == idx][plot_col].values
        if len(bin_data) > 0:
            boxplot_data.append(bin_data)
            boxplot_positions.append(idx)
    
    if len(boxplot_data) > 0:
        bp = ax.boxplot(boxplot_data, positions=boxplot_positions, widths=0.6,
                        patch_artist=False, showcaps=True, showfliers=False,
                        boxprops=dict(linewidth=1.0, color="black"),
                        whiskerprops=dict(linewidth=1.0, color="black"),
                        medianprops=dict(linewidth=1.0, color="black"),
                        capprops=dict(linewidth=0))
    
    # Highlight specific sample if requested
    if highlight_sample and highlight_sample != "":
        sample_df = chr_df[chr_df["sample"] == highlight_sample].sort_values("bin_idx")
        if len(sample_df) > 0:
            ax.plot(sample_df["bin_idx"], sample_df[plot_col], color="magenta", linewidth=2, label=" " + highlight_sample)
            ax.legend(loc="upper right", framealpha=0.9)
    
    # Format axes
    ax.set_xlabel(f"{chromosome} Position (Binned)", fontsize=12)
    if chromosome in ["chrX", "chrY"] or sex_subset is None:
        ax.set_ylabel("Normalized Depth", fontsize=12)
    else:
        ax.set_ylabel("Estimated Copy Number", fontsize=12)
    ax.set_xticks([])  # Remove bin tick labels as in R script
    
    # Add sample count text
    n_samples = chr_df["sample"].nunique()
    ax.text(0.02, 0.98, f"N={n_samples} samples",
           transform=ax.transAxes, va="top", ha="left", fontsize=10,
           bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    
    plt.tight_layout()
    
    # Save with appropriate filename
    filename = f"estimated_CN_per_bin.{suffix}.{chromosome}.png"
    raw_depth_dir = os.path.join(output_dir, "raw_depth_binned")
    os.makedirs(raw_depth_dir, exist_ok=True)
    output_path = os.path.join(raw_depth_dir, filename)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Created plot: {raw_depth_dir}/{filename}")


def generate_estimatePloidy_plots(df, bin_stats_df, sex_df, output_dir, highlight_sample=""):
    """
    Generate all plots that replicate estimatePloidy.R output.
    
    Args:
        df: DataFrame with chromosome-level statistics
        bin_stats_df: DataFrame with bin-level statistics (optional)
        sex_df: DataFrame with sex assignments (sample, sex, chrX_CN, chrY_CN)
        output_dir: Directory to save plots
        highlight_sample: Optional sample ID to highlight in plots
        
    Returns:
        None
    """
    print("\nGenerating estimatePloidy.R-style plots...")
    
    # 1. Sex assignment scatter plot
    plot_sex_assignments(sex_df, output_dir, highlight_sample)
    
    # 2. Per-contig boxplots - males
    print("\nGenerating per-contig boxplots for males...")
    plot_cn_per_contig_boxplot(df, output_dir, sex_subset="male", connect_samples=True, highlight_sample=highlight_sample)
    plot_cn_per_contig_boxplot(df, output_dir, sex_subset="male", connect_samples=False, highlight_sample=highlight_sample)
    
    # 3. Per-contig boxplots - females
    print("\nGenerating per-contig boxplots for females...")
    plot_cn_per_contig_boxplot(df, output_dir, sex_subset="female", connect_samples=True, highlight_sample=highlight_sample)
    plot_cn_per_contig_boxplot(df, output_dir, sex_subset="female", connect_samples=False, highlight_sample=highlight_sample)
    
    # 4. Per-contig boxplots - all samples
    print("\nGenerating per-contig boxplots for all samples...")
    plot_cn_per_contig_boxplot(df, output_dir, sex_subset="all", connect_samples=True, highlight_sample=highlight_sample)
    plot_cn_per_contig_boxplot(df, output_dir, sex_subset="all", connect_samples=False, highlight_sample=highlight_sample)
    
    # 5. Per-bin boxplots for each chromosome (if bin_stats_df provided)
    if bin_stats_df is not None and len(bin_stats_df) > 0:
        print("\nGenerating per-bin boxplots for each chromosome...")
        
        # Autosomes for all samples
        autosomes = [f"chr{i}" for i in range(1, 23)]
        for chrom in autosomes:
            if chrom in bin_stats_df["chr"].unique():
                plot_cn_per_bin_chromosome(bin_stats_df, output_dir, chrom, sex_subset=None, highlight_sample=highlight_sample)
        
        # Sex chromosomes for all samples
        if "chrX" in bin_stats_df["chr"].unique():
            plot_cn_per_bin_chromosome(bin_stats_df, output_dir, "chrX", sex_subset=None, highlight_sample=highlight_sample)
        
        if "chrY" in bin_stats_df["chr"].unique():
            plot_cn_per_bin_chromosome(bin_stats_df, output_dir, "chrY", sex_subset=None, highlight_sample=highlight_sample)
    
    print("\nestimaPloidy.R-style plots complete!")


def generate_aneuploidy_plots(df, bin_plot_data, training_loss_df, output_dir):
    """
    Generate all aneuploidy detection plots.
    
    Args:
        df: DataFrame with chromosome-level statistics
        bin_plot_data: DataFrame with bin-level data for plotting
        training_loss_df: DataFrame with training loss history
        output_dir: Directory to save plots
        
    Returns:
        None
    """
    print("\nGenerating aneuploidy detection plots...")
    
    # Plot training loss
    print("  Plotting training loss...")
    plot_training_loss(training_loss_df, output_dir)
    
    # Plot bin variance and bias
    print("  Plotting bin variance and bias...")
    plot_bin_variance_bias(bin_plot_data, output_dir)
    
    # Get all sample variance values
    all_sample_vars = bin_plot_data.groupby("sample")["sample_var"].first().values
    
    # Separate samples into aneuploid and normal
    aneuploid_samples = df[df["is_aneuploid"]]["sample"].unique()
    normal_samples = df[~df["sample"].isin(aneuploid_samples)]["sample"].unique()
    
    print(f"  Skipping per-sample plots (found {len(aneuploid_samples)} aneuploid and {len(normal_samples)} normal samples)")
    print("  Aneuploidy detection plots complete!")


def parse_args():
    """
    Parse command-line arguments.

    Args:
        None

    Returns:
        argparse.Namespace object with parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Aggregate and visualize aneuploidy detection results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "input",
        help="Input TSV file with aneuploidy detection results",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Output directory for plots and aneuploid TSV",
    )
    parser.add_argument(
        "-p",
        "--truth-json",
        required=False,
        default=None,
        help="JSON file with sample to aneuploidy type mappings (NORMAL, TRISOMY_13, TRISOMY_18, TRISOMY_21, KLINEFELTER, TRIPLE_X, TURNER, JACOBS, MULTIPLE, OTHER)",
    )
    parser.add_argument(
        "-e",
        "--exclusion-list",
        required=False,
        default=None,
        help="Optional text file with list of exclusion case IDs (one per line)",
    )
    parser.add_argument(
        "-b",
        "--bin-stats",
        required=False,
        default=None,
        help="Optional bin-level statistics TSV from aneuploidy_pyro.py (bin_stats.tsv.gz)",
    )
    parser.add_argument(
        "-t",
        "--training-loss",
        required=False,
        default=None,
        help="Optional training loss TSV from aneuploidy_pyro.py (training_loss.tsv)",
    )
    parser.add_argument(
        "--highlight-sample",
        required=False,
        default="",
        help="Optional sample ID to highlight in estimatePloidy.R-style plots",
    )
    parser.add_argument(
        "--skip-per-sample-plots",
        action="store_true",
        help="Skip generating individual plots for each sample (speeds up analysis for large cohorts)",
    )

    return parser.parse_args()


def main():
    """
    Main function to aggregate and visualize aneuploidy data.

    Args:
        None

    Returns:
        None
    """
    args = parse_args()

    # Create output directory
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")

    # Read TSV file
    print(f"Reading TSV file: {args.input}")
    df = pd.read_csv(args.input, sep="\t")

    # Display basic statistics
    print("\nData summary:")
    print(f"  Total samples: {df['sample'].nunique()}")
    print(f"  Total rows: {len(df)}")
    print(
        f"  Aneuploid rows: {df['is_aneuploid'].sum()} ({100 * df['is_aneuploid'].sum() / len(df):.2f}%)"
    )
    print(f"  Chromosomes: {sorted(df['chromosome'].unique())}")

    # Load truth and predict aneuploidy_types
    print("\nPredicting aneuploidy types...")
    if args.truth_json:
        truth_dict = load_truth_json(args.truth_json)
        print(f"  Loaded {len(truth_dict)} known truth cases")
    else:
        truth_dict = {}
        print("  No truth JSON provided - all samples will be compared against NORMAL")

    # Exclude cases if provided
    if args.exclusion_list:
        print(f"\nExcluding exclusion cases from {args.exclusion_list}...")
        exclusion_ids = set(load_exclusion_ids(args.exclusion_list))

        # Filter dataframe to exclude samples
        n_samples_before = df["sample"].nunique()
        df = df[~df["sample"].isin(exclusion_ids)]
        n_samples_after = df["sample"].nunique()
        print(f"  Excluded {n_samples_before - n_samples_after} samples from analysis")
        print(f"  Remaining samples: {n_samples_after}")

    print("\nAssigning sex and aneuploidy types...")
    aneuploidy_type_df = assign_sex_and_aneuploidy_types(df, truth_dict)

    if args.truth_json:
        print("Calculating prediction accuracy metrics...")
        calculate_metrics(aneuploidy_type_df, output_dir)
    else:
        print("Skipping metrics calculation (no truth JSON provided)")

    # Generate plots
    print("Generating histograms by chromosome type...")
    plot_histograms_by_chr_type(df, output_dir, args.highlight_sample)

    plot_aneuploid_histograms(df, output_dir, args.highlight_sample)

    # Generate median depth distribution plots
    generate_median_depth_distributions(df, output_dir, args.highlight_sample)
    
    # Generate estimatePloidy.R-style plots if bin stats is provided
    if args.bin_stats:
        print("\nLoading bin statistics for estimatePloidy.R-style plots...")
        bin_stats_df = pd.read_csv(args.bin_stats, sep="\t", compression="gzip")
        generate_estimatePloidy_plots(df, bin_stats_df, aneuploidy_type_df, output_dir, args.highlight_sample)
    else:
        print("\nGenerating estimatePloidy.R-style plots (chromosome-level only, no bin-level plots)...")
        generate_estimatePloidy_plots(df, None, aneuploidy_type_df, output_dir, args.highlight_sample)
    
    # Generate aneuploidy detection plots if bin stats is provided
    if args.bin_stats and args.training_loss:
        if args.skip_per_sample_plots:
            print("\nSkipping aneuploidy detection plots (--skip-per-sample-plots enabled)")
        else:
            print("\nLoading bin statistics and training loss for aneuploidy detection plots...")
            bin_stats_df = pd.read_csv(args.bin_stats, sep="\t", compression="gzip")
            training_loss_df = pd.read_csv(args.training_loss, sep="\t")
            generate_aneuploidy_plots(df, bin_stats_df, training_loss_df, output_dir)
    elif args.bin_stats or args.training_loss:
        print("\nWarning: Both --bin-stats and --training-loss are required for aneuploidy detection plots")
        print("Skipping aneuploidy detection plots...")
    else:
        print("\nSkipping aneuploidy detection plots (no bin statistics or training loss provided)")

    # Export aneuploid data
    aneuploid_output = os.path.join(output_dir, "aneuploid_samples.tsv")
    export_aneuploid_data(df, aneuploid_output)

    print("\nDone!")


if __name__ == "__main__":
    main()
