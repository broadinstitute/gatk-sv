import argparse
import os
import sys
import tempfile
from statistics import median
import gzip
import pybedtools
import pandas as pd
from pybedtools import BedTool


## part 1: (from line 1 to Line 94 of process_posthoc_cpx_depth_regenotyping.sh)
# reassign variant labels based on depth regenotyping in mod04b
def parse_arguments():
    """Parse command-line arguments."""

    # Create the parser
    parser = argparse.ArgumentParser(description="Reassign variant labels based on depth regenotyping in mod04b")

    # Positional arguments
    parser.add_argument("INVCF", help="Original input VCF prior to regenotyping")
    parser.add_argument("INTERVALS", help="BED file of genotyped intervals")
    parser.add_argument("GENOTYPES", help="Melted depth genotypes")
    parser.add_argument("FAMFILE", help=".fam file for all samples to be considered")
    parser.add_argument("OUTVCF", help="Full path to output VCF after relabeling")

    # Optional arguments
    parser.add_argument("-s", "--MINSIZE", type=int, default=1000,
                        help="Minimum size (in bp) of CNV interval to be considered during variant reclassification [default: 1000]")
    parser.add_argument("-d", "--MINDIFF", type=float, default=0.4,
                        help="Minimum difference in non-ref CNV genotype frequency between predicted carriers and noncarriers to consider a CNV interval to have adequate depth support [default: 0.4]")
    parser.add_argument("-D", "--MINSIZEiDEL", type=int, default=150,
                        help="Minimum insertion site size (in bp) to be considered for distinguishing insertion site deletions [default: 150 bp]")
    parser.add_argument("-T", "--MINdDUPTHRESH", type=int, default=1000000,
                        help="Minimum size (in bp) at which to prioritize an inverted dDUP classification over a dupINV or INVdup classification [default: 1000000 bp]")
    parser.add_argument("-R", "--RTABLE", default=None,
                        help="Path to table containing the final reclassification decision made per variant. [default: no table output]")
    parser.add_argument("-G", "--GTCOUNTSTABLE", default=None,
                        help="Path to table containing the raw genotype counts table per interval per variant. [default: no table output]")

    # Parse the arguments
    args = parser.parse_args()

    return args


## part 2: (from line 94 to line 172 of process_posthoc_cpx_depth_regenotyping.sh)

def validate_inputs(args):
    """Validate input arguments and files."""

    # Validate INVCF
    if not args.INVCF or not os.path.exists(args.INVCF):
        print("\nERROR: input VCF either empty or not found\n")
        sys.exit(0)
    if not args.INVCF.endswith('.gz'):
        print("\nERROR: input VCF must be bgzipped\n")
        sys.exit(0)

    # Validate INTERVALS
    if not args.INTERVALS or not os.path.exists(args.INTERVALS):
        print("\nERROR: intervals file either empty or not found\n")
        sys.exit(0)
    if not args.INTERVALS.endswith('.gz'):
        print("\nERROR: intervals file must be bgzipped\n")
        sys.exit(0)

    # Validate GENOTYPES
    if not args.GENOTYPES or not os.path.exists(args.GENOTYPES):
        print("\nERROR: melted genotypes file either empty or not found\n")
        sys.exit(0)
    if not args.GENOTYPES.endswith('.gz'):
        print("\nERROR: melted genotypes file must be bgzipped\n")
        sys.exit(0)

    # Validate FAMFILE
    if not args.FAMFILE or not os.path.exists(args.FAMFILE):
        print("\nERROR: input .fam file either empty or not found\n")
        sys.exit(0)

    # Validate OUTVCF
    if not args.OUTVCF:
        print("\nERROR: path to output VCF not specified\n")
        sys.exit(0)

    # Validate MINSIZE
    if args.MINSIZE is None:
        print("\nERROR: minimum CNV size not set\n")
        sys.exit(0)

    # Validate MINDIFF
    if args.MINDIFF is None:
        print("\nERROR: minimum CNV carrier freq difference not set\n")
        sys.exit(0)


def get_execution_directory():
    """Return the path to the execution directory."""
    return os.path.dirname(os.path.realpath(__file__))


def create_temp_directory():
    """Create and return a temporary directory."""
    return tempfile.mkdtemp()


## part 3: (from line 172 to line 223 of process_posthoc_cpx_depth_regenotyping.sh)
###PREPARE SAMPLE lists

def prepare_sample_lists(fam_file, output_dir):
    """Prepare sample lists based on the .fam file."""

    # Function to extract samples based on gender code
    def extract_samples_by_gender(lines, gender_code):
        return [line.split()[1] for line in lines if line[4] == gender_code]

    # Check if the file is gzipped and read accordingly
    if fam_file.endswith('.gz'):
        import gzip
        with gzip.open(fam_file, 'rt') as f:
            lines = [line for line in f if not line.startswith("#")]
    else:
        with open(fam_file, 'r') as f:
            lines = [line for line in f if not line.startswith("#")]

    # All samples
    with open(os.path.join(output_dir, 'all.samples.list'), 'w') as f:
        for sample in set(line.split()[1] for line in lines):
            f.write(sample + '\n')

    # Male samples
    with open(os.path.join(output_dir, 'male.samples.list'), 'w') as f:
        for sample in extract_samples_by_gender(lines, '1'):
            f.write(sample + '\n')

    # Female samples
    with open(os.path.join(output_dir, 'female.samples.list'), 'w') as f:
        for sample in extract_samples_by_gender(lines, '2'):
            f.write(sample + '\n')

###GATHER DATA PER VARIANT
def setup_genotype_counts_header(output_dir):
    """Set up header for genotype_counts_per_variant.bed."""
    header = "#chr\tstart\tend\tVID\tcarrier_del\tcarrier_wt\tcarrier_dup\tcontrol_del\tcontrol_wt\tcontrol_dup\tdiff_case_control_del_frac\tdiff_case_control_dup_frac"
    with open(os.path.join(output_dir, 'genotype_counts_per_variant.bed'), 'w') as f:
        f.write(header + '\n')

# del medCN # Not sure if the variable exists


## part 4: (from line 223 to line 312 of process_posthoc_cpx_depth_regenotyping.sh)

###Get list of samples & reference CN to consider, dependent on chr of call
# Summary: get_sample_lists_by_chromosome function processes sample lists based on the chromosome in question
# and decides which samples should be treated as carriers and which should be treated as controls.
# It also determines the default copy number (medCNdefault) based on the chromosome.
#ChrX: use only diploid females if possible, otherwise use haploid males

def get_sample_lists_by_chromosome(chr_name, samps, output_dir):
    """
    Decide which sample list to use based on the chromosome.

    Parameters:
    - chr_name (str): Name of the chromosome (e.g., "X", "chrY").
    - samps (str): Comma-separated list of samples.
    - output_dir (str): Directory where the sample lists are stored.

    Returns:
    - int: The default copy number for the given chromosome.
    """

    # Define paths for the carrier and control sample temporary files
    carrier_samples = os.path.join(output_dir, 'carrier_samples.tmp')
    control_samples = os.path.join(output_dir, 'control_samples.tmp')

    # Inner function to write a list of samples to a specified file
    def write_samples_to_file(samples, file_path):
        with open(file_path, 'w') as f:
            for sample in samples:
                f.write(sample + '\n')

    # If the chromosome is X or chrX
    if chr_name in ["X", "chrX"]:
        # Load the list of female samples
        female_samples = set(open(os.path.join(output_dir, 'female.samples.list')).read().split())
        # Find the intersection between the provided samples and the female samples
        matched_female_carriers = set(samps.split(',')) & female_samples
        if matched_female_carriers:
            # Write the matched female carriers to the carrier samples file
            write_samples_to_file(matched_female_carriers, carrier_samples)
            # Write the unmatched female samples to the control samples file
            write_samples_to_file(female_samples - matched_female_carriers, control_samples)
            # Default copy number for X chromosome in females
            medCNdefault = 2
        else:
            # If no matched female carriers, load the list of male samples
            male_samples = set(open(os.path.join(output_dir, 'male.samples.list')).read().split())
            # Find the intersection between the provided samples and the male samples
            matched_male_carriers = set(samps.split(',')) & male_samples
            # Write the matched male carriers and unmatched male samples to the respective files
            write_samples_to_file(matched_male_carriers, carrier_samples)
            write_samples_to_file(male_samples - matched_male_carriers, control_samples)
            # Default copy number for X chromosome in males
            medCNdefault = 1

    # If the chromosome is Y or chrY
    elif chr_name in ["Y", "chrY"]:
        # Load the list of male samples
        male_samples = set(open(os.path.join(output_dir, 'male.samples.list')).read().split())
        # Find the intersection between the provided samples and the male samples
        matched_male_carriers = set(samps.split(',')) & male_samples
        # Write the matched male carriers and unmatched male samples to the respective files
        write_samples_to_file(matched_male_carriers, carrier_samples)
        write_samples_to_file(male_samples - matched_male_carriers, control_samples)
        # Default copy number for Y chromosome
        medCNdefault = 1
    else:
        # For all other chromosomes
        # Load the list of all samples
        all_samples = set(open(os.path.join(output_dir, 'all.samples.list')).read().split())
        # Find the intersection between the provided samples and all samples
        matched_carriers = set(samps.split(',')) & all_samples
        # Write the matched carriers and unmatched samples to the respective files
        write_samples_to_file(matched_carriers, carrier_samples)
        write_samples_to_file(all_samples - matched_carriers, control_samples)
        # Default copy number for autosomal chromosomes
        medCNdefault = 2

    return medCNdefault

###Get genotype counts per variant
# Summary: get_genotype_counts_per_variant function processes the genotypes file and
# counts the number of carriers and controls for each variant.
# (gathers genotypes from a genotypes file for a specific genomic interval,
# and the method to do this depends on the size of the interval.)
def gather_genotypes_by_interval(chr_name, start, end, VID, genotypes, output_dir):
    """
    Gather genotypes corresponding to the interval of interest.

    Parameters:
    - chr_name (str): Name of the chromosome (e.g., "X", "chrY").
    - start (int): Start position of the interval.
    - end (int): End position of the interval.
    - VID (str): Variant identifier.
    - genotypes (str): Path to the genotypes file.
    - output_dir (str): Directory where the output will be stored.

    Returns:
    - str: Path to the temporary file containing intersected genotypes.
    """

    # Define path for the temporary file containing intersected genotypes
    intersected_genotypes_path = os.path.join(output_dir, 'intersected_genotypes.tmp.bed')

    # Convert genotypes to a BedTool object
    genotypes_bed = pybedtools.BedTool(genotypes)

    # Create a BedTool object for the interval of interest
    interval_bed = pybedtools.BedTool(f"{chr_name} {start} {end}", from_string=True)

    if end - start < 1000000:
        # Intersect the genotypes with the interval and filter by VID
        intersected_genotypes = genotypes_bed.intersect(interval_bed, wa=True, r=True, f=0.95)
        df = pybedtools.BedTool.to_dataframe(intersected_genotypes)
        df = df[df['name'] == VID]
    else:
        # Convert genotypes to a pandas DataFrame and filter by VID
        df = pybedtools.BedTool.to_dataframe(genotypes_bed)
        df = df[df['name'] == VID]

        # Compute coverage (this step will keep only the rows overlapping with the interval)
        df = df[(df['start'] <= end) & (df['end'] >= start)]
        df = df[['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart']]

    # Save the results to a file
    df.to_csv(intersected_genotypes_path, sep='\t', header=False, index=False)

    return intersected_genotypes_path


## part 5: (from line 312 to line 409 of process_posthoc_cpx_depth_regenotyping.sh)

def calculate_median_CN(control_samples_file, intersected_genotypes_file):
    """
    Calculate the median copy number across predicted non-carrier controls.

    Parameters:
    - control_samples_file (str): Path to the file containing the list of control samples.
    - intersected_genotypes_file (str): Path to the file containing intersected genotypes.

    Returns:
    - float or None: Median copy number value for the control samples, or None if there are no values.
    """

    # Read control samples from the file into a list
    with open(control_samples_file, 'r') as f:
        control_samples = [line.strip() for line in f]

    # Initialize an empty list to store copy number values from control samples
    CN_values = []

    # Read the intersected genotypes file line by line
    with open(intersected_genotypes_file, 'r') as f:
        for line in f:
            parts = line.split('\t')
            sample, CN = parts[3], parts[5]

            # If the sample from the current line is in the control samples list,
            # add its copy number value to the CN_values list
            if sample in control_samples:
                CN_values.append(int(CN))

    # Return the median of the copy number values, or None if the list is empty
    return median(CN_values) if CN_values else None


def count_genotypes_by_median(sample_file, intersected_genotypes_file, medCN):
    """Count genotypes lower than, equal to, or greater than the overall median."""
    with open(sample_file, 'r') as f:
        samples = [line.strip() for line in f]

    lower, equal, greater = 0, 0, 0
    with open(intersected_genotypes_file, 'r') as f:
        for line in f:
            parts = line.split('\t')
            sample, CN = parts[3], int(parts[5])
            if sample in samples:
                if CN < medCN:
                    lower += 1
                elif CN == medCN:
                    equal += 1
                else:
                    greater += 1
    return lower, equal, greater

def process_genotypes(carrier_samples_file, control_samples_file, intersected_genotypes_file, intervals_file, output_file, medCN):
    """
    Process genotypes to count and compute statistics based on median copy number (medCN).

    Parameters:
    - carrier_samples_file (str): Path to the file containing carrier samples.
    - control_samples_file (str): Path to the file containing control samples.
    - intersected_genotypes_file (str): Path to the intersected genotypes file.
    - intervals_file (str): Path to the intervals file.
    - output_file (str): Path to the output file.
    - medCN (float): Median copy number.

    Returns:
    None
    """

    # Load the intersected genotypes into a DataFrame
    df_genotypes = pd.read_csv(intersected_genotypes_file, sep='\t', header=None, names=['chr', 'start', 'end', 'sample', 'score', 'CN', 'thickStart'])

    # Load the carrier and control samples into lists
    with open(carrier_samples_file, 'r') as f:
        carrier_samples = [line.strip() for line in f]
    with open(control_samples_file, 'r') as f:
        control_samples = [line.strip() for line in f]

    # Count genotypes for carriers and non-carriers
    results = []
    for sample_list, label in [(carrier_samples, 'carrier'), (control_samples, 'control')]:
        if sample_list:
            sub_df = df_genotypes[df_genotypes['sample'].isin(sample_list)]
            less_than_medCN = len(sub_df[sub_df['CN'] < medCN])
            equal_to_medCN = len(sub_df[sub_df['CN'] == medCN])
            greater_than_medCN = len(sub_df[sub_df['CN'] > medCN])
            results.extend([less_than_medCN, equal_to_medCN, greater_than_medCN])
        else:
            results.extend([0, 0, 0])

    # Compute statistics
    ncase = results[2] + results[3] + results[4]
    nctrl = results[5] + results[6] + results[7]
    nctrl = 1 if nctrl == 0 else nctrl
    if ncase > 0:
        stats1 = (results[0] / ncase) - (results[5] / nctrl)
        stats2 = (results[2] / ncase) - (results[7] / nctrl)
    else:
        stats1, stats2 = 0, 0

    # Append results to the output file
    with open(output_file, 'a') as f:
        f.write("\t".join(map(str, results + [stats1, stats2])) + "\n")

# Further translation for the interval cleaning and processing will be done next.


def clean_intervals(output_dir, MINSIZE, INTERVALS):
    """
    Clean up intervals, process intervals per variant, and write results.

    Parameters:
    - output_dir (str): Directory where the output will be stored.
    - MINSIZE (int): Minimum size threshold for intervals.
    - INTERVALS (str): Path to the intervals file.

    Returns:
    None
    """
    header = "#chr\tstart\tend\tinterval_size\tinterval_type\tVID\tvariant_svtype\tvariant_cpxtype\tcarrier_del\tcarrier_wt\tcarrier_dup\tcontrol_del\tcontrol_wt\tcontrol_dup\tdiff_case_control_del_frac\tdiff_case_control_dup_frac\tCNV_assessment"

    cleaned_intervals_file = os.path.join(output_dir, 'genotype_counts_per_variant.cleaned_intervals.bed')

    # Read the genotype counts per variant into a DataFrame
    df_genotypes = pd.read_csv(cleaned_intervals_file, sep='\t', comment='#')

    # Extract unique mergedVIDs
    unique_VIDs = df_genotypes['VID'].unique()

    results = []

    for mergedVID in unique_VIDs:
        # Parse the mergedVID
        VID, cpxtype, svtype, *rest_intervals = mergedVID.split(';')

        if rest_intervals == ["NA"]:
            # If intervals are "NA", get them from the INTERVALS file
            df_intervals = pd.read_csv(INTERVALS, sep='\t', comment='#')
            intervals = df_intervals[df_intervals['VID'] == mergedVID]['interval'].tolist()
        else:
            intervals = rest_intervals

        for interval in intervals:
            interval_type, chr_name, start_end = interval.split('_')
            start, end = map(int, start_end.split('-'))

            # Filter rows from df_genotypes that overlap with the current interval
            overlapping_rows = df_genotypes[
                (df_genotypes['chr'] == chr_name) &
                (df_genotypes['start'] <= end) &
                (df_genotypes['end'] >= start)
                ]

            for _, row in overlapping_rows.iterrows():
                # Calculate the assessment based on the logic provided
                diff_case = row['diff_case_control_del_frac']
                diff_ctrl = row['diff_case_control_dup_frac']
                if diff_case > MINDIFF and diff_ctrl > MINDIFF:
                    assess = "DELDUP"
                elif diff_case > MINDIFF:
                    assess = "DEL"
                elif diff_ctrl > MINDIFF:
                    assess = "DUP"
                elif end - start < MINSIZE:
                    assess = "TOO_SMALL"
                else:
                    assess = "WT"

                results.append(
                    [chr_name, start, end, end - start, interval_type, VID, svtype, cpxtype] + row.values.tolist() + [
                        assess])

    # Convert the results list to a DataFrame
    df_results = pd.DataFrame(results)

    # Sort the DataFrame
    df_results = df_results.sort_values(
        by=[0, 1, 2, 4, 3])  # Sorting by columns: chr, start, end, interval_type, interval_size

    # Write the results to the cleaned_intervals_file
    df_results.to_csv(cleaned_intervals_file, sep='\t', header=False, index=False)

## part 6: (from line 409 to line 489 of process_posthoc_cpx_depth_regenotyping.sh)
# Convert Relevant Variants from VCF to BED:
#
# 1- Extract variants from the VCF that are not single-ended inversions and convert them to BED format using svtk vcf2bed.
# 2- Extract variants that are single-ended inversions and convert them to BED format.
# 3- Then, use modified coordinates for the inversion single-enders.

# Open and read the VCF file
def read_vcf(file_path):
    with gzip.open(file_path, 'rt') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        filepath_or_buffer = file_path,
        header = None,
        delimiter = '\t',
        comment = '#'
    )

# Read in the cleaned_intervals.bed file
# Convert results to a DataFrame and save to a text file
final_assessment = pd.DataFrame(results)
final_assessment.to_csv("${GTDIR}/final_variant_reassessment_table.txt", sep="\t", header=True, index=False)



def convert_variants_to_bed(INVCF, output_dir):
    """Convert relevant variants from VCF to BED format."""
    variants_to_reclassify_file = os.path.join(output_dir, 'variants_to_reclassify.vcf2bed.bed')
    cleaned_intervals_file = os.path.join(output_dir, 'genotype_counts_per_variant.cleaned_intervals.bed')

    header_lines = []
    relevant_variants = []

    with gzip.open(INVCF, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                header_lines.append(line.strip())
            else:
                fields = line.strip().split("\t")
                INFO_field = fields[7]
                if "SVTYPE=INV" not in INFO_field:
                    # Not an INV single ender, add to relevant_variants list
                    relevant_variants.append(line.strip())
                else:
                    # Handle INV single enders and modify their coordinates
                    chrom, start, _, _, _, _, _, info = fields
                    end = info.split("END=")[1].split(";")[0]
                    new_variant = f"{chrom}\t{start}\t{end}"
                    relevant_variants.append(new_variant)

    with open(variants_to_reclassify_file, 'w') as f:
        f.write("\n".join(header_lines + relevant_variants))

    return variants_to_reclassify_file

# Placeholder execution for the function
output_dir = "/path/to/output/directory"
INVCF = "/path/to/input/INVCF.vcf.gz"
convert_variants_to_bed(INVCF, output_dir)

# convert_variants_to_bed(INVCF, GTDIR)
# Convert to BED
# NOTE: The svtk vcf2bed conversion is more complex.
# For now, I'm keeping this step as a placeholder.
bed_data = "svtk vcf2bed conversion here"

# Save to BED file
bed_data.to_csv("${GTDIR}/variants_to_reclassify.vcf2bed.bed", sep="\t", header=False, index=False)


# Make Final Assessment for Each Variant:
#
# A header is printed for the final assessment table.
# The script retrieves the index of the SOURCE column from the BED file.
# Each structural variant (SV) is then evaluated. The classification of the SV is based on its type (svtype) and its predicted complex type (cpxtype).
# Depending on the SV's classification, it is processed in different ways, including:
# Complex variants (CPX) are assessed based on their predicted CNV intervals.
# Insertion variants (INS) are classified based on their complex type, e.g., translocational insertions, inverted dispersed duplications, etc.
# Single-ended inversions (BND) are checked for confirmation of palindromic inverted duplications.
# Other SV types are kept as they are with a reason of "IRRELEVANT_SV_TYPE".
# Initialize an empty list to store results
import os

def final_assessment_for_variant(output_dir):
    """Make a final assessment for each variant based on various criteria."""
    cleaned_intervals_file = os.path.join(output_dir, 'genotype_counts_per_variant.cleaned_intervals.bed')
    reassessment_table_file = os.path.join(output_dir, 'final_variant_reassessment_table.txt')

    # Print header
    with open(reassessment_table_file, 'w') as f:
        f.write(
            "#VID\tMODIFICATION\tREASON\tNEW_SVTYPE\tNEW_CPX_TYPE\tNEW_CPX_INTERVALS\tNEW_SVLEN\tNEW_SOURCE\tNEW_START\tNEW_END\n")

    # Read unique VID, svtype, and cpxtype from cleaned_intervals_file
    variant_data = {}
    with open(cleaned_intervals_file, 'r') as f:
        for line in f:
            if not line.startswith("#"):  # Assuming # are headers or comments
                parts = line.strip().split("\t")
                VID = parts[0]
                svtype = parts[1]  # Assuming svtype is the second column
                cpxtype = parts[2]  # Assuming cpxtype is the third column
                variant_data[VID] = (svtype, cpxtype)

    for VID, (svtype, cpxtype) in variant_data.items():
        MOD = "None"
        REASON = "None"
        if svtype == "CPX":
            # Logic to check CNV intervals and decide MOD and REASON
            # For this example, let's assume based on length
            NEW_SVLEN = int(cpxtype.split("_")[1])  # Assuming cpxtype has format "TypeX_LENGTH"
            if NEW_SVLEN > 1000:
                MOD = "Significant"
                REASON = "Length > 1000"
            else:
                MOD = "Minor"
                REASON = "Length <= 1000"

        NEW_SVTYPE = svtype
        NEW_CPX_TYPE = cpxtype
        NEW_CPX_INTERVALS = "INTERVALS_VALUE"  # Placeholder
        NEW_SOURCE = "SOURCE_VALUE"  # Placeholder
        NEW_START = "START_VALUE"  # Placeholder
        NEW_END = "END_VALUE"  # Placeholder

        # Write results to reassessment_table_file
        with open(reassessment_table_file, 'a') as f:
            f.write(f"{VID}\t{MOD}\t{REASON}\t{NEW_SVTYPE}\t{NEW_CPX_TYPE}\t{NEW_CPX_INTERVALS}\t{NEW_SVLEN}\t{NEW_SOURCE}\t{NEW_START}\t{NEW_END}\n")

# Placeholder execution for the


## part 7: (from line 489 to line 571 of process_posthoc_cpx_depth_regenotyping.sh)
# This function provides the structure for examining CNV for insertion.


def examine_cnv_for_insertion(VID, output_dir, MINSIZEiDEL):
    GTDIR = output_dir
    svtype = ""
    cpxtype = ""
    MOD = ""
    REASON = ""
    cpxintervals = ""
    SOURCE = ""
    SVLEN = 0

    def get_source_index(directory):
        with open(os.path.join(directory, "variants_to_reclassify.vcf2bed.bed"), 'r') as f:
            header = f.readline().strip().split('\t')
            return header.index('SOURCE') if 'SOURCE' in header else None

    SOURCEidx = get_source_index(GTDIR)

    import pandas as pd

    # Assuming a TSV file format
    def get_svtype_and_cpxtype_from_file(VID, directory):
        filepath = os.path.join(directory, "variant_data.tsv")
        df = pd.read_csv(filepath, sep='\t')

        # Filtering based on VID
        filtered_row = df[df['VID'] == VID]

        if not filtered_row.empty:
            svtype = filtered_row['svtype'].values[0]
            cpxtype = filtered_row['cpxtype'].values[0]
            return svtype, cpxtype
        else:
            return None, None

    # Inside the `examine_cnv_for_insertion` function
    svtype, cpxtype = get_svtype_and_cpxtype_from_file(VID, GTDIR)

    if svtype == "INS":
        if cpxtype in ["INS_A2B", "INS_B2A", "CTX_*INS_A2B", "CTX_*INS_B2A"]:
            bed = pybedtools.BedTool(os.path.join(GTDIR, "variants_to_reclassify.vcf2bed.bed"))
            sourcecoord_df = bed.filter(lambda b: b.name == VID).to_dataframe()
            sourcecoord = sourcecoord_df.iloc[0, SOURCEidx].replace("_", "\t").split("\t")[1] if not sourcecoord_df.empty else None

            if sourcecoord:
                cleaned_intervals_bed = pybedtools.BedTool(os.path.join(GTDIR, "genotype_counts_per_variant.cleaned_intervals.bed"))
                source_interval = pybedtools.create_interval_from_list(sourcecoord.replace(":", "\t").replace("-", "\t").split("\t"))
                sourceisdup = len(cleaned_intervals_bed.intersect(source_interval, wa=True, r=True, f=0.95).filter(lambda b: b.fields[-1] == "DUP"))

                sink_df = cleaned_intervals_bed.intersect(source_interval, wa=True, v=True, r=True, f=0.95).to_dataframe()
                sinkcoord = f"{sink_df.iloc[0, 0]}:{sink_df.iloc[0, 1]}-{sink_df.iloc[0, 2]}" if not sink_df.empty else None
                sinksize = int(sink_df.iloc[0, 3]) if not sink_df.empty else None

                if sinksize and sinksize >= MINSIZEiDEL:
                    sinkisdel = len(cleaned_intervals_bed.intersect(source_interval, wa=True, v=True, r=True, f=0.95).filter(lambda b: b.fields[-1] == "DEL"))
                else:
                    sinkisdel = 0

                if sourceisdup > 0 and sinkisdel > 0:
                    svtype = "CPX"
                    cpxtype = "dDUP_iDEL"
                    MOD = "RECLASSIFY"
                    REASON = "DISPERSED_DUPLICATION_W_INSERT_SITE_DEL"
                    cpxintervals = f"DUP_{sourcecoord},DEL_{sinkcoord}"
                    SOURCE = f"DUP_{sourcecoord}"
                    SVLEN = compute_svlen(cpxintervals)
                elif sourceisdup > 0 and sinkisdel == 0:
                    svtype = "CPX"
                    cpxtype = "dDUP"
                    MOD = "RECLASSIFY"
                    REASON = "DISPERSED_DUPLICATION"
                    cpxintervals = f"DUP_{sourcecoord}"
                    SOURCE = f"DUP_{sourcecoord}"
                elif sourceisdup == 0 and sinkisdel > 0:
                    svtype = "CPX"
                    cpxtype = "INS_iDEL"
                    MOD = "RECLASSIFY"
                    REASON = "INSERT_SITE_DEL"
                    cpxintervals = f"DEL_{sinkcoord}"
                    SOURCE = f"INS_{sourcecoord}"
                    SVLEN = compute_svlen(f"{cpxintervals},{sourcecoord}")
                else:
                    MOD = "KEEP"
                    REASON = "NON_DUPLICATED_INSERTION"
            else:
                MOD = "SKIP"
                REASON = "NO_SOURCE_INTERVAL_IN_CURRENT_SHARD"

    def compute_svlen(cpx_intervals):
        intervals = cpx_intervals.split(',')
        total_length = 0
        for interval in intervals:
            start, end = [int(x) for x in interval.split('_')[-1].split(':')[1].split('-')]
            total_length += end - start
        return total_length

    return {
        'svtype': svtype,
        'cpxtype': cpxtype,
        'MOD': MOD,
        'REASON': REASON,
        'cpxintervals': cpxintervals,
        'SOURCE': SOURCE,
        'SVLEN': SVLEN
    }

# Example call
result = examine_cnv_for_insertion("VID_example", "/path/to/output/dir", 150)
print(result)


## part 8: (from line 571 to line 677 of process_posthoc_cpx_depth_regenotyping.sh)


def examine_dup5_ins3(VID, output_dir, MINSIZEiDEL, MINdDUPTHRESH):
    GTDIR = output_dir
    SOURCEidx = get_source_index(GTDIR)
    svtype = ""
    cpxtype = ""
    MOD = ""
    REASON = ""
    cpxintervals = ""
    SOURCE = ""
    SVLEN = 0
    START = 0
    END = 0

    def get_source_index(directory):
        with open(os.path.join(directory, "variants_to_reclassify.vcf2bed.bed"), 'r') as f:
            header = f.readline().strip().split('\t')
            return header.index('SOURCE') if 'SOURCE' in header else None

    bed = pybedtools.BedTool(os.path.join(GTDIR, "variants_to_reclassify.vcf2bed.bed"))
    dupinterval_df = bed.filter(lambda b: b.name == VID).to_dataframe()
    dupinterval = dupinterval_df.iloc[0, SOURCEidx].replace("_", "\t").split("\t")[
        1] if not dupinterval_df.empty else None

    if dupinterval:
        dupsize = int(dupinterval.split(":")[1].split("-")[1]) - int(dupinterval.split(":")[1].split("-")[0])

        cleaned_intervals_bed = pybedtools.BedTool(
            os.path.join(GTDIR, "genotype_counts_per_variant.cleaned_intervals.bed"))
        dup_interval = pybedtools.create_interval_from_list(
            dupinterval.replace(":", "\t").replace("-", "\t").split("\t"))
        dupconfirmed = len(cleaned_intervals_bed.intersect(dup_interval, wa=True, r=True, f=0.95).filter(
            lambda b: b.fields[-1] != "WT"))

        if dupconfirmed > 0:
            sink_df = bed.filter(lambda b: b.name == VID).to_dataframe()
            sinkinterval = f"{sink_df.iloc[0, 0]}:{sink_df.iloc[0, 1]}-{sink_df.iloc[0, 2]}"
            sinksize = int(sinkinterval.split(":")[1].split("-")[1]) - int(sinkinterval.split(":")[1].split("-")[0])

            sink_interval_bed = pybedtools.create_interval_from_list(
                sinkinterval.replace(":", "\t").replace("-", "\t").split("\t"))
            sinkisdel = len(cleaned_intervals_bed.intersect(sink_interval_bed, wa=True, r=True, f=0.95).filter(
                lambda b: b.fields[-1] == "DEL"))

            inv_start = min(int(dupinterval.split(":")[1].split("-")[0]), int(sinkinterval.split(":")[1].split("-")[0]))
            inv_end = max(int(dupinterval.split(":")[1].split("-")[1]), int(sinkinterval.split(":")[1].split("-")[1]))
            invinterval = f"{dupinterval.split(':')[0]}:{inv_start}-{inv_end}"
            invsize = inv_end - inv_start

            if invsize < MINdDUPTHRESH:
                if sinkisdel > 0:
                    invinterval_bed = pybedtools.create_interval_from_list(
                        invinterval.replace(":", "\t").replace("-", "\t").split("\t"))
                    invinterval = invinterval_bed.subtract(sink_interval_bed).to_dataframe()
                    invinterval = f"{invinterval.iloc[0, 0]}:{invinterval.iloc[0, 1]}-{invinterval.iloc[0, 2]}"

                    svtype = "CPX"
                    cpxtype = "dupINVdel"
                    MOD = "RECLASSIFY"
                    REASON = "DUP_FLANKED_INVERSION_WITH_DEL"
                    cpxintervals = f"DUP_{dupinterval},INV_{invinterval},DEL_{sinkinterval}"
                    SVLEN = invsize
                    START = min(
                        [int(coord.split(":")[1].split("-")[0]) for coord in [dupinterval, invinterval, sinkinterval]])
                    END = max(
                        [int(coord.split(":")[1].split("-")[1]) for coord in [dupinterval, invinterval, sinkinterval]])
                else:
                    svtype = "CPX"
                    cpxtype = "dupINV"
                    MOD = "RECLASSIFY"
                    REASON = "DUP_FLANKED_INVERSION"
                    cpxintervals = f"DUP_{dupinterval},INV_{invinterval}"
                    SVLEN = invsize
                    START = min([int(coord.split(":")[1].split("-")[0]) for coord in [dupinterval, invinterval]])
                    END = max([int(coord.split(":")[1].split("-")[1]) for coord in [dupinterval, invinterval]])
            else:
                if sinkisdel > 0:
                    svtype = "CPX"
                    cpxtype = "dDUP_iDEL"
                    MOD = "RECLASSIFY"
                    REASON = "INVERTED_DISPERSED_DUPLICATION_WITH_DELETION"
                    cpxintervals = f"DUP_{dupinterval},INV_{dupinterval},DEL_{sinkinterval}"
                    SOURCE = f"DUP_{dupinterval}"
                    SVLEN = dupsize + sinksize
                    START = min(
                        [int(coord.split(":")[1].split("-")[0]) for coord in [dupinterval, invinterval, sinkinterval]])
                    END = max(
                        [int(coord.split(":")[1].split("-")[1]) for coord in [dupinterval, invinterval, sinkinterval]])
                else:
                    svtype = "CPX"
                    cpxtype = "dDUP"
                    MOD = "RECLASSIFY"
                    REASON = "INVERTED_DISPERSED_DUPLICATION"
                    cpxintervals = f"DUP_{dupinterval},INV_{dupinterval}"
                    SOURCE = f"DUP_{dupinterval}"
                    SVLEN = dupsize
                    START = int(dupinterval.split(":")[1].split("-")[0])
                    END = int(dupinterval.split(":")[1].split("-")[1])
        else:
            MOD = "UNRESOLVED"
            REASON = "PREDICTED_DUP_INTERVAL_FAILED_GT"

    return {
        'svtype': svtype,
        'cpxtype': cpxtype,
        'MOD': MOD,
        'REASON': REASON,
        'cpxintervals': cpxintervals,
        'SOURCE': SOURCE,
        'SVLEN': SVLEN,
        'START': START,
        'END': END
    }


# Example call
result = examine_dup5_ins3("VID_example", "/path/to/output/dir", 150, 1000)
print(result)

## part 9: (from line 677 to line 840 of process_posthoc_cpx_depth_regenotyping.sh)

def examine_cnv_evidence(VID, output_dir, MINSIZEiDEL, MINdDUPTHRESH, MINSIZE):
    # Directory with generated files
    GTDIR = output_dir
    # Variables initialization
    SOURCEidx = get_source_index(GTDIR)
    svtype = ""
    cpxtype = ""
    MOD = ""
    REASON = ""
    cpxintervals = ""
    SOURCE = ""
    SVLEN = 0
    START = 0
    END = 0

    # Nested function to get the index of 'SOURCE' from the header of a BED file
    def get_source_index(directory):
        with open(os.path.join(directory, "variants_to_reclassify.vcf2bed.bed"), 'r') as f:
            header = f.readline().strip().split('\t')
            return header.index('SOURCE') if 'SOURCE' in header else None

    def get_dup_values(VID, GTDIR, SOURCEidx):
        # Load the BED file into a pandas DataFrame
        bed_file_path = os.path.join(GTDIR, "variants_to_reclassify.vcf2bed.bed")
        df = pd.read_csv(bed_file_path, sep='\t', header=None)

        # Extract dupinterval
        dupinterval_row = df[df[0] == VID].iloc[0]
        dupinterval = dupinterval_row[SOURCEidx].split('_')[1]

        # Calculate dupsize
        start, end = map(int, dupinterval.split('-'))
        dupsize = end - start

        # Calculate dupconfirmed
        cleaned_intervals_path = os.path.join(GTDIR, "genotype_counts_per_variant.cleaned_intervals.bed")
        cleaned_df = pd.read_csv(cleaned_intervals_path, sep='\t', header=None)
        matching_intervals = cleaned_df[cleaned_df[0] == VID]
        confirmed_intervals = matching_intervals[~matching_intervals.iloc[:, -1].str.contains("WT")]
        dupconfirmed = len(confirmed_intervals)

        return dupinterval, dupsize, dupconfirmed

    def get_sink_values(VID, GTDIR):
        # Load the BED file into a pandas DataFrame
        bed_file_path = os.path.join(GTDIR, "variants_to_reclassify.vcf2bed.bed")
        df = pd.read_csv(bed_file_path, sep='\t', header=None)

        # Extract sinkinterval
        sinkinterval_row = df[df[0] == VID].iloc[0]
        sinkinterval = f"{sinkinterval_row[0]}:{sinkinterval_row[1]}-{sinkinterval_row[2]}"

        # Calculate sinksize
        start, end = sinkinterval_row[1], sinkinterval_row[2]
        sinksize = end - start

        # Calculate sinkisdel
        cleaned_intervals_path = os.path.join(GTDIR, "genotype_counts_per_variant.cleaned_intervals.bed")
        cleaned_df = pd.read_csv(cleaned_intervals_path, sep='\t', header=None)
        matching_intervals = cleaned_df[cleaned_df[0] == VID]
        del_intervals = matching_intervals[matching_intervals.iloc[:, -1].str.contains("DEL")]
        sinkisdel = len(del_intervals)

        return sinkinterval, sinksize, sinkisdel

    # Fetch the values for dupinterval, dupsize, and dupconfirmed
    dupinterval, dupsize, dupconfirmed = get_dup_values(VID, GTDIR, SOURCEidx)

    # Fetch the values for sinkinterval, sinksize, and sinkisdel
    sinkinterval, sinksize, sinkisdel = get_sink_values(VID, GTDIR)

    # Rest of the logic from the bash script
    if dupconfirmed > 0:
        if invsize < MINdDUPTHRESH:
            if sinkisdel > 0:
                # If sink is deleted, classify as delINVdup
                svtype = "CPX"
                cpxtype = "delINVdup"
                MOD = "RECLASSIFY"
                REASON = "DUP_FLANKED_INVERSION_WITH_DEL"
                cpxintervals = f"DEL_{sinkinterval},INV_{invinterval},DUP_{dupinterval}"
                SVLEN = invsize
                # START and END calculations (still need further logic for accurate computation)
            else:
                # If sink is not clearly deleted (or too small), classify as dupINV (no del)
                svtype = "CPX"
                cpxtype = "INVdup"
                MOD = "RECLASSIFY"
                REASON = "DUP_FLANKED_INVERSION"
                cpxintervals = f"INV_{invinterval},DUP_{dupinterval}"
                SVLEN = invsize
                # START and END calculations (still need further logic for accurate computation)
        else:
            if sinkisdel > 0:
                # If sink is deleted, classify as dDUP_iDEL
                svtype = "CPX"
                cpxtype = "dDUP_iDEL"
                MOD = "RECLASSIFY"
                REASON = "INVERTED_DISPERSED_DUPLICATION_WITH_DELETION"
                cpxintervals = f"DUP_{dupinterval},INV_{dupinterval},DEL_{sinkinterval}"
                SOURCE = f"DUP_{dupinterval}"
                SVLEN = dupsize + sinksize
            else:
                # If sink is not clearly deleted (or too small), classify as dDUP (no iDEL)
                svtype = "CPX"
                cpxtype = "dDUP"
                MOD = "RECLASSIFY"
                REASON = "INVERTED_DISPERSED_DUPLICATION"
                cpxintervals = f"INV_{dupinterval},DUP_{dupinterval}"
                SOURCE = f"DUP_{dupinterval}"
                SVLEN = dupsize
    else:
        MOD = "UNRESOLVED"
        REASON = "PREDICTED_DUP_INTERVAL_FAILED_GT"

    return {
        'svtype': svtype,
        'cpxtype': cpxtype,
        'MOD': MOD,
        'REASON': REASON,
        'cpxintervals': cpxintervals,
        'SOURCE': SOURCE,
        'SVLEN': SVLEN,
        'START': START,
        'END': END
    }


# Example call
result = examine_cnv_evidence("VID_example", "/path/to/output/dir", 150, 1000, 100)
print(result)



## part 10: (from line 840 to line 925 of process_posthoc_cpx_depth_regenotyping.sh)

def salvage_inverted_duplications(variant_ids, output_dir, RTABLE="", GTCOUNTSTABLE=""):
    # Directory with generated files
    GTDIR = output_dir
    # List to store results
    results = []

    for VID in variant_ids:
        # Variables initialization
        svtype = ""
        cpxtype = ""
        MOD = ""
        REASON = ""
        cpxintervals = ""
        SOURCE = ""
        SVLEN = 0
        START = 0
        END = 0

        # Extracting `dupinterval` from file
        with open(os.path.join(GTDIR, 'variants_to_reclassify.vcf2bed.bed'), 'r') as f:
            for line in f:
                if VID in line:
                    parts = line.strip().split()
                    dupinterval = f"{parts[0]}:{parts[1]}-{parts[2]}"
                    break

        # Calculating `dupsize`
        interval_parts = dupinterval.split(":")[1].split("-")
        SVLEN = int(interval_parts[1]) - int(interval_parts[0])

        # Extracting `dupconfirmed` from file
        with open(os.path.join(GTDIR, 'genotype_counts_per_variant.cleaned_intervals.bed'), 'r') as f:
            dupconfirmed = sum(1 for line in f if VID in line and "DUP" in line)

        # Logic based on `cpxtype`
        if cpxtype in ["INVERSION_SINGLE_ENDER_--", "INVERSION_SINGLE_ENDER_++"]:
            if dupconfirmed > 0:
                svtype = "CPX"
                MOD = "RECLASSIFY"
                REASON = "PALINDROMIC_INVERTED_DUPLICATION"
                cpxintervals = f"DUP_{dupinterval},INV_{dupinterval}"
                START, END = sorted(interval_parts)
                cpxtype = "piDUP_FR" if cpxtype == "INVERSION_SINGLE_ENDER_--" else "piDUP_RF"
            else:
                MOD = "KEEP"
                REASON = "DID_NOT_RESOLVE_AS_piDUP"
        else:
            MOD = "KEEP"
            REASON = "IRRELEVANT_SV_TYPE"

        # Finalize unset variables
        cpxintervals = cpxintervals if cpxintervals else "."
        SVLEN = SVLEN if SVLEN else "."
        SOURCE = SOURCE if SOURCE else "."
        START = START if START else "."
        END = END if END else "."

        # Append result for this variant
        results.append((VID, MOD, REASON, svtype, cpxtype, cpxintervals, SVLEN, SOURCE, START, END))

    # Write results to file
    with open(os.path.join(GTDIR, 'final_variant_reassessment_table.txt'), 'w') as f:
        for result in results:
            f.write('\t'.join(map(str, result)) + '\n')

    # Copy final reclassification table, if optioned
    if RTABLE:
        os.system(f"cp {os.path.join(GTDIR, 'final_variant_reassessment_table.txt')} {RTABLE}")

    # Copy genotyping counts, if optioned
    if GTCOUNTSTABLE:
        os.system(f"cp {os.path.join(GTDIR, 'genotype_counts_per_variant.cleaned_intervals.bed')} {GTCOUNTSTABLE}")

# Example call
variant_ids_example = ["VID1", "VID2"]  # Sample list of variant IDs
salvage_inverted_duplications(variant_ids_example, "/path/to/output/dir")

## part 11: (from line 925 to line 994 of process_posthoc_cpx_depth_regenotyping.sh)


def generate_final_vcf(input_vcf, reassessment_table, output_dir):
    GTDIR = output_dir
    reassessed_variants = {}  # To store reassessed variants for later processing

    # Prepare the CPX type descriptions
    cpx_types = [
        "##CPX_TYPE_delINV=\"Complex inversion with 5' flanking deletion.\"",
        "##CPX_TYPE_INVdel=\"Complex inversion with 3' flanking deletion.\"",
        #... (other CPX type descriptions) ...
    ]

    with gzip.open(input_vcf, 'rt') as invcf, open(os.path.join(GTDIR, 'cleaned_output.vcf'), 'w') as outvcf:
        # Extract header
        for line in invcf:
            if line.startswith("##"):
                outvcf.write(line)
            else:
                break  # Stop when header section is over

        # Add CPX type descriptions
        for description in cpx_types:
            outvcf.write(description + '\n')

        # Add sample line to VCF header (the first non-"##" line)
        outvcf.write(line)

        # Get variant IDs for reassessment
        variant_ids = set()
        with open(reassessment_table, 'r') as f:
            for line in f:
                if not line.startswith("#"):
                    variant_ids.add(line.split('\t')[0])

        # Write variants not in reassessment table to output VCF
        for line in invcf:
            variant_id = line.split('\t')[2]  # Get variant ID from the line
            if variant_id not in variant_ids:
                outvcf.write(line)
            else:
                reassessed_variants[variant_id] = line  # Store the variant for reassessment

    # Process reassessed variants
    with open(os.path.join(GTDIR, 'final_variant_reassessment_table.txt'), 'r') as reassessment_table:
        for line in reassessment_table:
            if line.startswith("#"):
                continue  # Skip header lines

            VID, MOD, REASON, svtype, cpxtype, cpxintervals, SVLEN, SOURCE, START, END = line.strip().split('\t')

            variant_line = reassessed_variants.get(VID)  # Get the original variant line
            if not variant_line:
                continue  # If the variant is not found, skip to the next

            fields = variant_line.strip().split('\t')

            if MOD in ["SKIP", "KEEP"]:
                with open(os.path.join(GTDIR, 'variants_to_be_reassessed.vcf'), 'a') as varvcf:
                    varvcf.write(variant_line)

            elif MOD == "UNRESOLVED":
                fields[6] = "UNRESOLVED"
                fields[7] += ";UNRESOLVED_TYPE=POSTHOC_RD_GT_REJECTION"
                with open(os.path.join(GTDIR, 'variants_to_be_reassessed.vcf'), 'a') as varvcf:
                    varvcf.write('\t'.join(fields) + '\n')

            elif MOD == "RECLASSIFY":
                if SVLEN == "." or START == "." or END == ".":
                    info_fields = {kv.split('=')[0]: kv.split('=')[1] for kv in fields[7].split(';') if '=' in kv}
                    SVLEN = info_fields.get('SVLEN', SVLEN)
                    START = fields[1]
                    END = info_fields.get('END', END)

                fields[1] = START
                new_info = ";".join([f"SVTYPE={svtype}", f"SVLEN={SVLEN}", f"END={END}"])
                fields[7] = new_info

                with open(os.path.join(GTDIR, 'variants_to_be_reassessed.vcf'), 'a') as varvcf:
                    varvcf.write('\t'.join(fields) + '\n')

    # To handle merging 'variants_to_be_reassessed.vcf' back into 'cleaned_output.vcf':
    with open(os.path.join(GTDIR, 'cleaned_output.vcf'), 'a') as outvcf:
        with open(os.path.join(GTDIR, 'variants_to_be_reassessed.vcf'), 'r') as varvcf:
            for line in varvcf:
                outvcf.write(line)

## part 12: (from line 925 to the end of of process_posthoc_cpx_depth_regenotyping.sh)

def refine_info_field(INFO, svtype, cpxtype, cpxintervals, SVLEN, SOURCE, END):
    # Use the logic from the Bash script to modify the INFO field
    INFO = INFO.replace(f"END=[^;]*;", f"END={END};")
    INFO = INFO.replace(f";END=[^;]*;", f";END={END};")
    INFO = INFO.replace(f";END=[^;]*$", f";END={END}")

    INFO = INFO.replace(f"SVTYPE=[^;]*;", f"SVTYPE={svtype};")
    INFO = INFO.replace(f";SVTYPE=[^;]*;", f";SVTYPE={svtype};")
    INFO = INFO.replace(f";SVTYPE=[^;]*$", f";SVTYPE={svtype}")

    INFO = INFO.replace(f"SVLEN=[^;]*;", f"SVLEN={SVLEN};")
    INFO = INFO.replace(f";SVLEN=[^;]*;", f";SVLEN={SVLEN};")
    INFO = INFO.replace(f";SVLEN=[^;]*$", f";SVLEN={SVLEN}")

    # Removing or modifying CPX_TYPE as needed
    if svtype == "CPX":
        if "CPX_TYPE=" not in INFO:
            INFO += f";CPX_TYPE={cpxtype}"
        else:
            INFO = INFO.replace(f";CPX_TYPE=[^;]*;", f";CPX_TYPE={cpxtype};")
            INFO = INFO.replace(f";CPX_TYPE=[^;]*$", f";CPX_TYPE={cpxtype}")
    else:
        INFO = INFO.replace(f";CPX_TYPE=[^;]*;", ";")
        INFO = INFO.replace(f";CPX_TYPE=[^;]*$", "")

    # Adding or modifying SOURCE as needed
    if SOURCE == ".":
        INFO = INFO.replace(f";SOURCE=[^;]*;", ";")
        INFO = INFO.replace(f";SOURCE=[^;]*$", "")
    else:
        if "SOURCE=" in INFO:
            INFO = INFO.replace(f";SOURCE=[^;]*;", f";SOURCE={SOURCE};")
            INFO = INFO.replace(f";SOURCE=[^;]*$", f";SOURCE={SOURCE}")
        else:
            INFO += f";SOURCE={SOURCE}"

    # Adding or modifying CPX_INTERVALS as needed
    if cpxintervals == ".":
        INFO = INFO.replace(f";CPX_INTERVALS=[^;]*;", ";")
        INFO = INFO.replace(f";CPX_INTERVALS=[^;]*$", "")
    else:
        if "CPX_INTERVALS=" in INFO:
            INFO = INFO.replace(f";CPX_INTERVALS=[^;]*;", f";CPX_INTERVALS={cpxintervals};")
            INFO = INFO.replace(f";CPX_INTERVALS=[^;]*$", f";CPX_INTERVALS={cpxintervals}")
        else:
            INFO += f";CPX_INTERVALS={cpxintervals}"

    return INFO


def refine_vcf_records(reassessment_table, variants_to_be_reassessed, cleaned_output_vcf, bin_dir, outvcf):
    with open(reassessment_table, 'r') as reassess, open(variants_to_be_reassessed, 'r') as var_vcf, open(
            cleaned_output_vcf, 'a') as out_vcf:
        vcf_lines = var_vcf.readlines()

        for line in reassess:
            if not line.startswith("#"):
                VID, MOD, REASON, svtype, cpxtype, cpxintervals, SVLEN, SOURCE, START, END = line.strip().split('\t')

                variant_record = next((line for line in vcf_lines if VID in line), None)
                if variant_record:
                    fields = variant_record.strip().split('\t')
                    INFO = refine_info_field(fields[7], svtype, cpxtype, cpxintervals, SVLEN, SOURCE, END)

                    # Print the modified record to the output VCF
                    out_vcf.write("\t".join(
                        [fields[0], START, fields[2], fields[3], "<{}>".format(svtype), fields[5], fields[6],
                         INFO] + fields[8:]) + "\n")

    # Sort and compress the final output VCF
    os.system(
        "{} {} | vcf-sort | bgzip -c > {}".format(os.path.join(bin_dir, 'rm_cpx_info.py'), cleaned_output_vcf, outvcf))

    # Clean up
    os.rmdir(GTDIR)

# Note: GTDIR variable is not provided in the code. It should be defined somewhere before calling os.rmdir(GTDIR).
# Note: cleaned_output_vcf variable is not provided in the code. It should be defined somewhere before calling os.system(...).
############################## end of part 12 of regenotyping pipeline
# This would be the main entry point if you run the Python script
if __name__ == "__main__":
    args = parse_arguments()

    # section 2
    validate_inputs(args)

    BIN = get_execution_directory()

    GTDIR = create_temp_directory()

    # Further processing can continue from here using the args object.
    invcf = args.INVCF
    INTERVALS = args.INTERVALS
    GENOTYPES = args.GENOTYPES
    FAMFILE = args.FAMFILE

    outvcf = args.OUTVCF
    #print("Input VCF: {}".format(invcf))
    #print("Output VCF: {}".format(outvcf))

    prepare_sample_lists(args.FAMFILE, GTDIR)

    setup_genotype_counts_header(GTDIR)
    # section 4:
    prepare_sample_lists(args.FAMFILE, GTDIR)

    setup_genotype_counts_header(GTDIR)

    # Assuming chr, start, end, VID, and samps are extracted from the loop reading variants
    # These variables should be set accordingly
    chr = "X"  # Placeholder
    start = 0  # Placeholder
    end = 1000  # Placeholder
    VID = "VID123"  # Placeholder
    samps = "sample1,sample2"  # Placeholder

    medCNdefault = get_sample_lists_by_chromosome(chr, samps, GTDIR)

    intersected_genotypes = gather_genotypes_by_interval(chr, start, end, VID, args.GENOTYPES, GTDIR)

