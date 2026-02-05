#!/bin/python

import argparse
import json
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.metrics import confusion_matrix, classification_report


def read_file_list(file_list_path):
    """Read a list of file paths from a text file."""
    with open(file_list_path, "r") as f:
        file_paths = [line.strip() for line in f if line.strip()]
    return file_paths


def concatenate_tsvs(file_paths):
    """Read and concatenate multiple TSV files into a single DataFrame."""
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
    """Classify chromosome as autosomal, chrX, or chrY."""
    if chrom == "chrX":
        return "chrX"
    elif chrom == "chrY":
        return "chrY"
    else:
        return "Autosomal"


def format_column_name(col):
    """Format column name for display."""
    return col.replace("_", " ").title()


def save_and_close_plot(output_dir, filename):
    """Save the current plot and close it."""
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, filename), dpi=300)
    plt.close()
    print(f"  Created histogram: {filename}")


def plot_single_histogram(data, x_col, title, xlabel, output_dir, filename, bins=50):
    """Create and save a single histogram plot."""
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
    data, x_col, hue_col, title, xlabel, output_dir, filename, palette, bins=50
):
    """Create and save a histogram plot with hue coloring."""
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
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Frequency")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    save_and_close_plot(output_dir, filename)


def plot_histograms_by_chr_type(df, output_dir):
    """
    Create histograms of all numeric columns (except specified exclusions),
    colored by chromosome type (autosomal, chrX, chrY).
    """
    # Add chromosome type column
    df = df.copy()
    df["chr_type"] = df["chromosome"].apply(get_chromosome_type)

    # Columns to exclude from plotting
    exclude_cols = ["sample", "chromosome", "chr_type"]

    # Get numeric columns
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    plot_cols = [col for col in numeric_cols if col not in exclude_cols]

    # Handle sample_var separately - convert to sample_stdev and plot by chr_type
    if "sample_var" in plot_cols:
        plot_cols.remove("sample_var")
        df["sample_stdev"] = np.sqrt(df["sample_var"])
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
        )


def plot_aneuploid_histograms(df, output_dir):
    """
    Create histograms for aneuploid samples only.
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
    plt.savefig(os.path.join(output_dir, "hist_aneuploid_chromosome.png"), dpi=300)
    plt.close()
    print("  Created histogram: hist_aneuploid_chromosome.png")

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
        )

    # Handle sample_stdev for aneuploid samples
    if "sample_var" in aneuploid_df.columns:
        aneuploid_df["sample_stdev"] = np.sqrt(aneuploid_df["sample_var"])
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


def predict_aneuploidy_type_from_aneuploidies(sample_df):
    """
    Predict aneuploidy_type based on aneuploidy calls for a sample.

    Args:
        sample_df: DataFrame with all chromosome rows for a single sample

    Returns:
        Predicted aneuploidy_type string
    """
    # Get all aneuploid chromosomes and their copy numbers
    aneuploid_calls = sample_df[sample_df["is_aneuploid"]]

    if len(aneuploid_calls) == 0:
        # No aneuploidies detected - determine sex based on chrX/chrY
        chrX_cn = (
            sample_df[sample_df["chromosome"] == "chrX"]["copy_number"].values[0]
            if "chrX" in sample_df["chromosome"].values
            else 2
        )
        chrY_cn = (
            sample_df[sample_df["chromosome"] == "chrY"]["copy_number"].values[0]
            if "chrY" in sample_df["chromosome"].values
            else 0
        )

        if chrX_cn == 2 and chrY_cn == 0:
            return "FEMALE"
        elif chrX_cn == 1 and chrY_cn == 1:
            return "MALE"
        else:
            return "UNKNOWN"

    if len(aneuploid_calls) > 1:
        return "MULTIPLE"

    # Single aneuploidy - determine aneuploidy_type
    chrom = aneuploid_calls.iloc[0]["chromosome"]
    cn = aneuploid_calls.iloc[0]["copy_number"]

    # Get chrX and chrY copy numbers for sex chromosome aneuploidies
    chrX_cn = (
        sample_df[sample_df["chromosome"] == "chrX"]["copy_number"].values[0]
        if "chrX" in sample_df["chromosome"].values
        else 2
    )
    chrY_cn = (
        sample_df[sample_df["chromosome"] == "chrY"]["copy_number"].values[0]
        if "chrY" in sample_df["chromosome"].values
        else 0
    )

    # Autosomal trisomies
    if chrom == "chr13" and cn == 3:
        return "TRISOMY_13"
    elif chrom == "chr18" and cn == 3:
        return "TRISOMY_18"
    elif chrom == "chr21" and cn == 3:
        return "TRISOMY_21"

    # Sex chromosome aneuploidies
    elif chrX_cn == 2 and chrY_cn == 1:
        return "KLINEFELTER"  # XXY
    elif chrX_cn == 3 and chrY_cn == 0:
        return "TRIPLE_X"  # XXX
    elif chrX_cn == 1 and chrY_cn == 0:
        return "TURNER"  # X0
    elif chrX_cn == 1 and chrY_cn == 2:
        return "JACOBS"  # XYY

    return "OTHER"


def assign_aneuploidy_types(df, truth_dict):
    """
    Assign true and predicted aneuploidy_types to each sample.

    Args:
        df: Combined dataframe with all samples and chromosomes
        truth_dict: Dictionary mapping sample IDs to true aneuploidy_types (empty dict if not provided)

    Returns:
        DataFrame with columns: sample, true_aneuploidy_type, predicted_aneuploidy_type, score
    """
    sample_aneuploidy_types = []

    # Group by sample for efficient iteration
    for sample_id, sample_df in df.groupby("sample"):
        # Get true aneuploidy_type from truth.json
        true_aneuploidy_type = truth_dict.get(sample_id, "NORMAL")

        # Create lookup dictionaries for faster access
        chrom_to_cn = dict(zip(sample_df["chromosome"], sample_df["copy_number"]))
        chrom_to_aneuploid = dict(
            zip(sample_df["chromosome"], sample_df["is_aneuploid"])
        )

        chrX_cn = chrom_to_cn.get("chrX", 2)
        chrY_cn = chrom_to_cn.get("chrY", 0)

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
        List of exclusion IDs
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

    Note:
        Converts sample_var to inferred_sample_std and removes sample_var and chr_type columns
    """
    aneuploid_df = df[df["is_aneuploid"]].copy()
    aneuploid_df["inferred_sample_std"] = np.sqrt(aneuploid_df["sample_var"])
    aneuploid_df = aneuploid_df.drop(columns=["sample_var"], errors="ignore")
    aneuploid_df = aneuploid_df.drop(columns=["chr_type"], errors="ignore")
    aneuploid_df.to_csv(output_path, sep="\t", index=False)
    print(f"\nExported {len(aneuploid_df)} aneuploid rows to {output_path}")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Aggregate and visualize aneuploidy detection results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input file containing list of TSV file paths (one per line)",
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

    return parser.parse_args()


def main():
    """Main function to aggregate and visualize aneuploidy data."""
    args = parse_args()

    # Create output directory
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")

    # Read file list and concatenate TSVs
    file_paths = read_file_list(args.input)
    df = concatenate_tsvs(file_paths)

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

    print("\nAssigning aneuploidy types...")
    aneuploidy_type_df = assign_aneuploidy_types(df, truth_dict)

    if args.truth_json:
        print("Calculating prediction accuracy metrics...")
        calculate_metrics(aneuploidy_type_df, output_dir)
    else:
        print("Skipping metrics calculation (no truth JSON provided)")

    # Generate plots
    print("Generating histograms by chromosome type...")
    plot_histograms_by_chr_type(df, output_dir)

    plot_aneuploid_histograms(df, output_dir)

    # Export aneuploid data
    aneuploid_output = os.path.join(output_dir, "aneuploid_samples.tsv")
    export_aneuploid_data(df, aneuploid_output)

    print("\nDone!")


if __name__ == "__main__":
    main()
