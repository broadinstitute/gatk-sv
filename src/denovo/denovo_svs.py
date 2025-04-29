"""
SV de novo filtering script
"""

import argparse
import yaml
import numpy as np
import pandas as pd
import pybedtools
import collections
import time
from subprocess import Popen, PIPE


pd.options.mode.chained_assignment = None  # default='warn'


def verbose_print(msg, verbose):
    if verbose == "True":
        print(msg)


def get_count(row, items):
    return sum(samp in items for samp in row['samples'].split(','))


def get_parents_frequency(row, items):
    return sum(s in items for s in row['samples'].split(',')) / len(items)


def get_family_count(row, ped):
    sample = row['sample']
    father = ped[(ped['IndividualID'] == sample)]['FatherID'].values[0]
    mother = ped[(ped['IndividualID'] == sample)]['MotherID'].values[0]
    parents = [father, mother]
    return (
        sum(s in parents for s in row['samples'].split(','))
    )


def variant_info(row, field, vcf):
    row_name = row['name']
    samp = row['sample']

    filt_pos = vcf[(vcf['ID'] == row_name)]['FORMAT'].str.split(':').tolist()

    if field in filt_pos[0]:
        idx = filt_pos[0].index(field)
        if samp in vcf.columns:
            samp_info = vcf[(vcf['ID'] == row_name)][samp].str.split(':').tolist()
            samp_info_field = samp_info[0][idx]
        else:
            samp_info_field = 'NA'
    else:
        samp_info_field = 'NA'

    return samp_info_field


def add_family(row, ped):
    sample = row['sample']
    fam = ped[(ped['IndividualID'] == sample)]['FamID'].values
    chr = row['chrom']
    if fam.size != 0:
        return '_'.join([str(fam[0]), chr])


def get_info_parent(row, ped, vcf, parent, field):
    row_name = row['name']
    sample = row['sample']
    parent_id = ped[(ped['IndividualID'] == sample)][parent].values

    filt_pos = vcf[(vcf['ID'] == row_name)]['FORMAT'].str.split(':').tolist()

    if field in filt_pos[0]:
        idx = filt_pos[0].index(field)
        parent_info = vcf[(vcf['ID'] == row_name)][parent_id[0]].str.split(':').tolist()
        parent_info_field = parent_info[0][idx]
    else:
        return 'NA'

    return parent_info_field


def get_parents(proband, ped):
    mother = ped[(ped['IndividualID'] == proband)]['MotherID'].values
    father = ped[(ped['IndividualID'] == proband)]['FatherID'].values
    parent_list = [mother, father]
    return parent_list


def get_family_id(row, ped):
    samp = row['sample']
    family_id = ped[(ped['IndividualID'] == samp)]['FamID'].values
    return family_id[0]


def get_bincov_matrix_url(sample, sample_batches, batch_bincov):
    batch = sample_batches.loc[sample_batches['sample'] == sample]['batch'].iloc[0]
    bincov_gs = batch_bincov.loc[batch_bincov['batch'] == batch]['bincov'].iloc[0]
    return bincov_gs


def write_to_filter_file(filename, header, original_df, filtered_df):
    filename.write(header)
    filename.write("\n")
    original = original_df['name_famid']
    filtered = filtered_df['name_famid']
    to_write = list(set(original) - set(filtered))
    filename.write(str(to_write))
    filename.write("\n")


def write_to_size_file(filename, header, df):
    filename.write(header)
    filename.write(str(len(df)))
    filename.write("\n")


def min_gq(row, gq_min):
    paternal_srgq = row['paternal_srgq']
    maternal_srgq = row['maternal_srgq']
    paternal_pegq = row['paternal_pegq']
    maternal_pegq = row['maternal_pegq']
    paternal_gq = row['paternal_gq']
    maternal_gq = row['maternal_gq']
    minimum = min([paternal_srgq, maternal_srgq, paternal_gq, maternal_gq, paternal_pegq, maternal_pegq])
    if minimum != '.':
        if int(min([paternal_srgq, maternal_srgq, paternal_gq, maternal_gq, paternal_pegq, maternal_pegq])) <= gq_min:
            return 'Remove'
        else:
            return 'Keep'
    else:
        return 'Keep'


def convert_to_bedtool(df, cols_to_keep=None, sort=True):
    if cols_to_keep is None:
        string_obj = df.to_string(header=False, index=False)
        bt_obj = pybedtools.BedTool(string_obj, from_string=True)
        if sort:
            string_obj = df.to_string(header=False, index=False)
            bt_obj = pybedtools.BedTool(string_obj, from_string=True).sort()
    else:
        string_obj = df[cols_to_keep].to_string(header=False, index=False)
        bt_obj = pybedtools.BedTool(string_obj, from_string=True).sort()
    return bt_obj


def tabix_query(filename, chrom, start, end):
    process = Popen(['tabix', '-h', filename, f'{chrom}:{start}-{end}'], stdout=PIPE)
    cov = []
    for line in process.stdout:
        cov.append(line.decode("utf-8").strip().split())
    cov_matrix = pd.DataFrame(cov)
    return cov_matrix


def get_per_sample_matrix(row, ped, sample_batches, batch_bincov, family_member):
    svtype = row['svtype']
    start = int(row['start'])
    end = int(row['end'])
    if svtype == 'INS':
        start = start - 100
        end = end + 100
    proband = row['sample']
    proband_matrix = get_bincov_matrix_url(proband, sample_batches, batch_bincov)
    mom = get_parents(proband, ped)[0][0]
    mom_matrix = get_bincov_matrix_url(mom, sample_batches, batch_bincov)
    dad = get_parents(proband, ped)[1][0]
    dad_matrix = get_bincov_matrix_url(dad, sample_batches, batch_bincov)
    if family_member == 'proband':
        return [proband, proband_matrix]
    if family_member == 'mother':
        return [mom, mom_matrix]
    if family_member == 'father':
        return [dad, dad_matrix]


def find_coverage(row, ped, sample_batches, batch_bincov, family_member):
    chrom = row['chrom']
    start = int(row['start'])
    end = int(row['end'])
    svtype = row['svtype']
    if svtype == 'INS':
        start = start - 100
        end = end + 100
    if svtype == 'DEL' or svtype == 'DUP' or svtype == 'INS':
        if family_member == 'mother':
            matrix = get_per_sample_matrix(row, ped, sample_batches, batch_bincov, family_member='mother')[1]
        if family_member == 'father':
            matrix = get_per_sample_matrix(row, ped, sample_batches, batch_bincov, family_member='father')[1]
        cov_matrix = tabix_query(matrix, chrom, start, end)

        header = cov_matrix.iloc[0].to_list()  # grab the first row for the header
        new_header = []
        for x in header:
            new_header.append(x)
        cov_matrix = cov_matrix[1:]  # take the data less the header row
        cov_matrix.columns = new_header
        if family_member == 'mother':
            indv = get_per_sample_matrix(row, ped, sample_batches, batch_bincov, family_member='mother')[0]
        if family_member == 'father':
            indv = get_per_sample_matrix(row, ped, sample_batches, batch_bincov, family_member='father')[0]
        sample_cov_matrix = cov_matrix[indv].tolist()
        sample_cov_matrix_decoded = []
        for x in sample_cov_matrix:
            sample_cov_matrix_decoded.append(x)
        return sample_cov_matrix_decoded
    else:
        return 'NA'


def get_median_coverage(mom_matrix, dad_matrix, coverage_cutoff):
    mom_coverage = mom_matrix
    dad_coverage = dad_matrix

    if mom_coverage == 'NA' and dad_coverage == 'NA':
        return 'Keep'

    if mom_coverage != 'NA' and dad_coverage != 'NA':
        coverage_d = {
            'Mother': mom_coverage,
            'Father': dad_coverage
        }

        coverage_df = pd.DataFrame(coverage_d)
        '''
        if coverage_df.empty:
            print('Coverage df is empty!')
            exit()
        '''
        median = coverage_df.median()
        median_list = median.tolist()
        median_filtered = [x for x in median_list if x < coverage_cutoff]
        if len(median_filtered) > 0:
            return 'Remove'
        else:
            return 'Keep'
    else:
        return 'Remove'  # TODO: check with Alba what should this return.


def get_cnv_intersection_depth(bed, raw, overlap):
    intersect = bed.coverage(raw).to_dataframe(disable_auto_names=True, header=None)
    if len(intersect) != 0:
        names_overlap = intersect[intersect[10] > overlap][6].to_list()
    else:
        names_overlap = ['']
    return names_overlap


def get_cnv_intersection_other(bed, raw, overlap):
    overlap = bed.intersect(raw, wo=True, f=overlap, r=True).to_dataframe(disable_auto_names=True, header=None)
    if len(overlap) != 0:
        names_overlap = overlap[6].to_list()
    else:
        names_overlap = ['']
    return names_overlap


def get_insertion_intersection(bed, raw):
    overlap = bed.closest(raw).to_dataframe(disable_auto_names=True, header=None)
    if len(overlap) != 0:
        overlap['is_close'] = abs(overlap[8] - overlap[1]) < 100
        if sum(bool(x) for x in overlap['is_close']) > 0:
            names_overlap = overlap[(overlap['is_close'])][6].to_list()
        else:
            names_overlap = ['']
    return names_overlap


def main():
    """
    Parse arguments vcf_metrics
    """
    parser = argparse.ArgumentParser(description='Parse arguments')
    parser.add_argument('--bed', dest='bed', help='Input BED file')
    parser.add_argument('--ped', dest='ped', help='Ped file')
    parser.add_argument('--vcf', dest='vcf', help='VCF file')
    parser.add_argument('--disorder', dest='disorder', help='Genomic disorder regions')
    parser.add_argument('--out', dest='out', help='Output file with all variants')
    parser.add_argument('--out_de_novo', dest='out_de_novo', help='Output file with only de novo variants')
    parser.add_argument('--raw_proband', dest='raw_proband', help='Directory with raw SV calls - output from m04')
    parser.add_argument('--raw_parents', dest='raw_parents', help='Directory with raw SV calls - output from m04')
    parser.add_argument('--raw_depth_proband', dest='raw_depth_proband', help='Directory with raw SV depth calls - output from m04')
    parser.add_argument('--raw_depth_parents', dest='raw_depth_parents', help='Directory with raw depth SV calls - output from m04')
    parser.add_argument('--config', dest='config', help='Config file')
    parser.add_argument('--exclude_regions', dest='exclude_regions', help='File containing regions with known somatic mutations')
    parser.add_argument('--coverage', dest='coverage', help='File with batch in first column respective coverage file in second column')
    parser.add_argument('--sample_batches', dest='sample_batches', help='File with samples in first column and their respective batch in second column')
    parser.add_argument('--verbose', dest='verbose', help='Verbosity')
    args = parser.parse_args()

    bed_file = args.bed
    ped_file = args.ped
    vcf_file = args.vcf
    disorder_file = args.disorder
    out_file = args.out
    de_novo_out_file = args.out_de_novo
    raw_file_proband = args.raw_proband
    raw_file_parent = args.raw_parents
    raw_file_depth_proband = args.raw_depth_proband
    raw_file_depth_parent = args.raw_depth_parents
    verbose = args.verbose
    config_file = args.config
    exclude_regions = args.exclude_regions
    coverage = args.coverage
    batches = args.sample_batches

    with open(config_file, "r") as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    large_cnv_size = int(config['large_cnv_size'])
    gnomad_col = config['gnomad_col']
    alt_gnomad_col = config['alt_gnomad_col']
    gnomad_af = float(config['gnomad_AF'])
    parents_af = float(config['parents_AF'])
    large_raw_overlap = float(config['large_raw_overlap'])
    small_raw_overlap = float(config['small_raw_overlap'])
    cohort_af = float(config['cohort_AF'])
    coverage_cutoff = float(config['coverage_cutoff'])
    depth_only_size = float(config['depth_only_size'])
    parents_overlap = float(config['parents_overlap'])
    gq_min = float(config['gq_min'])

    # Read files
    verbose_print('Reading Input Files', verbose)
    bed = pd.read_csv(bed_file, sep='\t').replace(np.nan, '', regex=True)
    bed = bed[(bed['samples'] != "")]
    bed.rename(columns={'#chrom': 'chrom'}, inplace=True)
    vcf = pd.read_csv(vcf_file, sep='\t')
    ped = pd.read_csv(ped_file, sep='\t')
    disorder = pd.read_csv(disorder_file, sep='\t', header=None)
    raw_bed_colnames = ['ID', 'start', 'end', 'svtype', 'sample']
    raw_bed_child = pd.read_csv(raw_file_proband, sep='\t', names=raw_bed_colnames, header=None).replace(np.nan, '', regex=True)
    raw_bed_parent = pd.read_csv(raw_file_parent, sep='\t', names=raw_bed_colnames, header=None).replace(np.nan, '', regex=True)
    raw_bed_depth_child = pd.read_csv(raw_file_depth_proband, sep='\t', names=raw_bed_colnames, header=None).replace(np.nan, '', regex=True)
    raw_bed_depth_parent = pd.read_csv(raw_file_depth_parent, sep='\t', names=raw_bed_colnames, header=None).replace(np.nan, '', regex=True)
    exclude_regions = pd.read_csv(exclude_regions, sep='\t').replace(np.nan, '', regex=True)
    bincov_colnames = ['batch', 'bincov', 'index']
    sample_batches_colnames = ['sample', 'batch']
    bincov = pd.read_csv(coverage, sep='\t', names=bincov_colnames, header=None).replace(np.nan, '', regex=True)
    sample_batches = pd.read_csv(batches, sep='\t', names=sample_batches_colnames, header=None).replace(np.nan, '', regex=True)

    #########################
    # REFORMAT AND ANNOTATE #
    #########################
    # Exit out of bed file is empty
    if bed.size == 0:
        bed.to_csv(path_or_buf=out_file, mode='a', index=False, sep='\t', header=True)
        bed.to_csv(path_or_buf=de_novo_out_file, mode='a', index=False, sep='\t', header=True)
        exit()
    # Get parents and children ids
    verbose_print('Getting parents/children/affected/unaffected IDs', verbose)
    start = time.time()
    list_of_trios = [item for item, count in collections.Counter(ped["FamID"]).items() if count == 3]
    families = ped[(ped['FamID'].isin(list_of_trios)) &
                   (ped['FatherID'] != "0") &
                   (ped['MotherID'] != "0")]['FamID'].values
    trios = ped[(ped['FamID'].isin(families))]
    parents = trios[(trios['FatherID'] == "0") & (trios['MotherID'] == "0")]['IndividualID'].values
    children = trios[(trios['FatherID'] != "0") & (trios['MotherID'] != "0")]['IndividualID'].values
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Get counts at cohort level and number of parents/children
    verbose_print('Getting counts', verbose)
    start = time.time()
    bed['num_children'] = bed.apply(lambda r: get_count(r, children), axis=1)
    bed['num_parents'] = bed.apply(lambda r: get_count(r, parents), axis=1)
    bed['AF_parents'] = bed.apply(lambda r: get_parents_frequency(r, parents), axis=1)
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Remove mCNVs, BNDs and SVs in sex chromosomes
    start = time.time()
    verbose_print('Remove BND and mCNV', verbose)
    bed = bed[(~bed['svtype'].isin(['BND', 'CNV']))]
    # bed = bed[(~bed['svtype'].isin(['BND', 'CNV']) & (~bed['chrom'].isin(["chrY", "chrX"])))]
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Flag if small or large CNV based on large_cnv_size cutoff and flag for removal based on size for depth only calls
    verbose_print('Flagging calls depending on size', verbose)
    start = time.time()
    bed['is_large_cnv'] = (bed['SVLEN'] >= large_cnv_size) & ((bed['svtype'] == 'DEL') | (bed['svtype'] == 'DUP'))
    bed['is_small_cnv'] = (bed['SVLEN'] < large_cnv_size) & ((bed['svtype'] == 'DEL') | (bed['svtype'] == 'DUP'))
    bed['is_depth_only'] = (bed['EVIDENCE'] == "RD")
    bed['is_depth_only_small'] = (bed['svtype'] == "DUP") & (bed['ALGORITHMS'] == "depth") & (bed['SVLEN'] <= depth_only_size)
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Split into one row per sample
    verbose_print('Split into one row per sample', verbose)
    start = time.time()
    bed_split = bed.assign(sample=bed.samples.str.split(",")).explode('sample')
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Sepparate variants in children and parents
    verbose_print('Sepparate variants in children and parents', verbose)
    start = time.time()
    bed_child = bed_split[bed_split['sample'].isin(children)]
    bed_parents = bed_split[bed_split['sample'].isin(parents)]
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Flag variants in genomic disorder (GD) regions
    verbose_print('Flagging variants in GD', verbose)
    start = time.time()
    bed_child['in_gd'] = bed_child['name'].isin(disorder[0])
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Annotate family information for filtering
    bed_child['family_id'] = bed_child.apply(lambda r: get_family_id(r, ped), axis=1)
    bed_child['name_famid'] = bed_child['name'] + "_" + bed_child['family_id'].astype(str).str.strip("[]")
    bed_parents['family_id'] = bed_parents.apply(lambda r: get_family_id(r, ped), axis=1)
    bed_parents['name_famid'] = bed_parents['name'] + "_" + bed_parents['family_id'].astype(str).str.strip("[]")

    # Filter out by frequency - AF gnomad < 0.01 OR inGD
    verbose_print('Filtering by frequency', verbose)
    start = time.time()
    bed_child["AF"] = pd.to_numeric(bed_child["AF"])
    try:
        bed_child[gnomad_col] = pd.to_numeric(bed_child[gnomad_col])
        remove_freq = bed_child[~((((bed_child[gnomad_col] <= gnomad_af) & (bed_child['AF'] <= cohort_af)) |
                                   ((bed_child[gnomad_col].isnull()) & (bed_child['AF'] <= cohort_af)) |
                                   (bed_child['in_gd'])))]['name_famid'].to_list()
        bed_child['is_de_novo'] = pd.Series(True, index=bed_child.index).mask(bed_child['name_famid'].isin(remove_freq), False)
        bed_child['filter_flag'] = pd.Series('de_novo', index=bed_child.index).mask(bed_child['name_famid'].isin(remove_freq), 'AF')

    except KeyError:
        bed_child[alt_gnomad_col] = pd.to_numeric(bed_child[alt_gnomad_col])
        remove_freq = bed_child[~((((bed_child[alt_gnomad_col] <= gnomad_af) & (bed_child['AF'] <= cohort_af)) |
                                   ((bed_child[alt_gnomad_col].isnull()) & (bed_child['AF'] <= cohort_af)) |
                                   (bed_child['in_gd'])))]['name_famid'].to_list()
        bed_child['is_de_novo'] = pd.Series(True, index=bed_child.index).mask(bed_child['name_famid'].isin(remove_freq), False)
        bed_child['filter_flag'] = pd.Series('de_novo', index=bed_child.index).mask(bed_child['name_famid'].isin(remove_freq), 'AF')
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Get counts within family and remove if SV in parents
    verbose_print('Keep variants in children only', verbose)
    start = time.time()
    try:
        parents_sc = int(config['parents_SC'])
    except KeyError:
        parents_sc = len(parents)
    if len(bed_child.index) > 0:  # Otherwise, if nothing in bed_child will get an error
        bed_child['num_parents_family'] = bed_child.apply(lambda r: get_family_count(r, ped), axis=1)

        remove_child = bed_child[~((bed_child['num_parents_family'] == 0) &  # there are no parents of the affected child in the samples column
                                   (bed_child['num_children'] >= 1) &  # for each line there is at least one child in the samples column
                                   (bed_child['AF_parents'] <= parents_af) &  # for each line, the frequency of parents in the samples column must be < 0.05
                                   (bed_child['num_parents'] <= parents_sc))]['name_famid'].to_list()  # for each line, the number of parents must be <= parents_SC or the total number of parents in the ped

        bed_child.loc[bed_child['name_famid'].isin(remove_child) & bed_child['is_de_novo'], 'filter_flag'] = 'in_parent'
        bed_child.loc[bed_child['name_famid'].isin(remove_child) & bed_child['is_de_novo'], 'is_de_novo'] = False
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Extract info from the VCF file
    verbose_print('Appending FILTER information', verbose)
    start = time.time()
    if len(bed_child.index) > 0:
        bed_child['GT'] = bed_child.apply(lambda r: variant_info(r, 'GT', vcf), axis=1)
        bed_child['EV'] = bed_child.apply(lambda r: variant_info(r, 'EV', vcf), axis=1)
        bed_child['GQ'] = bed_child.apply(lambda r: variant_info(r, 'GQ', vcf), axis=1)
        bed_child['RD_CN'] = bed_child.apply(lambda r: variant_info(r, 'RD_CN', vcf), axis=1)
        bed_child['RD_GQ'] = bed_child.apply(lambda r: variant_info(r, 'RD_GQ', vcf), axis=1)
        bed_child['PE_GQ'] = bed_child.apply(lambda r: variant_info(r, 'PE_GQ', vcf), axis=1)
        bed_child['PE_GT'] = bed_child.apply(lambda r: variant_info(r, 'PE_GT', vcf), axis=1)
        bed_child['SR_GQ'] = bed_child.apply(lambda r: variant_info(r, 'SR_GQ', vcf), axis=1)
        bed_child['SR_GT'] = bed_child.apply(lambda r: variant_info(r, 'SR_GT', vcf), axis=1)
    else:
        bed_child['GT'] = 'NA'
        bed_child['EV'] = 'NA'
        bed_child['GQ'] = 'NA'
        bed_child['RD_CN'] = 'NA'
        bed_child['RD_GQ'] = 'NA'
        bed_child['PE_GQ'] = 'NA'
        bed_child['PE_GT'] = 'NA'
        bed_child['SR_GQ'] = 'NA'
        bed_child['SR_GT'] = 'NA'
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Remove WHAM only and GT = 1
    verbose_print('Remove wham only and GT=1 calls', verbose)
    start = time.time()
    remove_wham = bed_child[(bed_child['ALGORITHMS'] == "wham") & (bed_child['GQ'] == '1')]
    bed_child.loc[bed_child['name_famid'].isin(remove_wham) & bed_child['is_de_novo'], 'filter_flag'] = 'wham_only'
    bed_child.loc[bed_child['name_famid'].isin(remove_wham), 'is_de_novo'] = False
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Annotate parental information
    if len(bed_child.index) > 0:
        bed_child['paternal_gq'] = bed_child.apply(lambda r: get_info_parent(r, ped, vcf, 'FatherID', 'GQ'), axis=1)
        bed_child['maternal_gq'] = bed_child.apply(lambda r: get_info_parent(r, ped, vcf, 'MotherID', 'GQ'), axis=1)
        bed_child['paternal_rdcn'] = bed_child.apply(lambda r: get_info_parent(r, ped, vcf, 'FatherID', 'RD_CN'), axis=1)
        bed_child['maternal_rdcn'] = bed_child.apply(lambda r: get_info_parent(r, ped, vcf, 'MotherID', 'RD_CN'), axis=1)
        bed_child['paternal_pegq'] = bed_child.apply(lambda r: get_info_parent(r, ped, vcf, 'FatherID', 'PE_GQ'), axis=1)
        bed_child['maternal_pegq'] = bed_child.apply(lambda r: get_info_parent(r, ped, vcf, 'MotherID', 'PE_GQ'), axis=1)
        bed_child['paternal_srgq'] = bed_child.apply(lambda r: get_info_parent(r, ped, vcf, 'FatherID', 'SR_GQ'), axis=1)
        bed_child['maternal_srgq'] = bed_child.apply(lambda r: get_info_parent(r, ped, vcf, 'MotherID', 'SR_GQ'), axis=1)
    else:
        bed_child['paternal_gq'] = 'NA'
        bed_child['maternal_gq'] = 'NA'
        bed_child['paternal_rdcn'] = 'NA'
        bed_child['maternal_rdcn'] = 'NA'
        bed_child['paternal_pegq'] = 'NA'
        bed_child['maternal_pegq'] = 'NA'
        bed_child['paternal_srgq'] = 'NA'
        bed_child['maternal_srgq'] = 'NA'

    # LARGE CNV: Check for false negative in parents: check depth in parents, independently of the calls
    verbose_print('Large CNVs check', verbose)
    start = time.time()
    # 1. Strip out if same CN in parents and proband if not in chrX
    remove_large_cnv = bed_child.loc[~(((bed_child['SVLEN'] <= 5000) | ((bed_child['RD_CN'] != bed_child['maternal_rdcn']) & (bed_child['RD_CN'] != bed_child['paternal_rdcn']))) | (bed_child['chrom'] == 'chrX'))]
    bed_child.loc[bed_child['name_famid'].isin(remove_large_cnv) & bed_child['is_de_novo'], 'filter_flag'] = 'same_cn_as_parents'
    bed_child.loc[bed_child['name_famid'].isin(remove_large_cnv) & bed_child['is_de_novo'], 'is_de_novo'] = False
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # 2. Check if call in parents with bedtools coverage (332)
    verbose_print('CNV present in parents check', verbose)
    start = time.time()
    cols_keep = ['family_chrom', 'start', 'end', 'name', 'svtype', 'sample', 'name_famid']
    bed_child_large = bed_child[(bed_child['is_large_cnv'])]
    bed_parents_large = bed_parents[(bed_parents['is_large_cnv'])]
    if (len(bed_child_large.index) > 0) & (len(bed_parents_large.index) > 0):
        bed_child_large['family_chrom'] = bed_child_large.apply(lambda r: add_family(r, ped), axis=1)
        bed_child_large = bed_child_large[cols_keep].to_string(header=False, index=False)
        bed_child_large = pybedtools.BedTool(bed_child_large, from_string=True).sort()

        bed_parents_large['family_chrom'] = bed_parents_large.apply(lambda r: add_family(r, ped), axis=1)
        bed_parents_large = bed_parents_large[cols_keep].to_string(header=False, index=False)
        bed_parents_large = pybedtools.BedTool(bed_parents_large, from_string=True).sort()

        bed_overlap = bed_child_large.coverage(bed_parents_large).to_dataframe(disable_auto_names=True, header=None)
        names_overlap = bed_overlap[(bed_overlap[10] >= parents_overlap)][6].to_list()
    else:
        names_overlap = ['']
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Add if overlap to bed child table
    bed_child['overlap_parent'] = (bed_child['name_famid'].isin(names_overlap))

    # Small calls:
    # If RD,SR and < large_cnv_size, treat RD,SR as SR
    verbose_print('Small CNVs check', verbose)
    start = time.time()
    bed_child['EVIDENCE_FIX'] = bed_child['EVIDENCE']
    bed_child[(bed_child['SVLEN'] <= large_cnv_size) &
              (bed_child['EVIDENCE'] == "RD,SR") &
              ((bed_child['svtype'] == 'DEL') | (bed_child['svtype'] == 'DUP'))]['EVIDENCE_FIX'] = "SR"

    # Add if variant contains RD evidence
    bed_child['contains_RD'] = bed_child.EVIDENCE_FIX.str.contains('RD')
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    ###################
    # RAW FILES CHECK #
    ###################
    # Keep all CPT and CTX - they are rare and the type is different than in the raw files
    # Check only DEL,DUP and INS for raw evidence
    verbose_print('Checking raw files', verbose)
    # Remove depth raw calls > 1MB
    raw_bed_depth_parent['SVLEN'] = raw_bed_depth_parent['end'] - raw_bed_depth_parent['start']
    raw_bed_ref_depth_parent_subset = raw_bed_depth_parent[(raw_bed_depth_parent['SVLEN'] < 1000000)]

    # Reformat raw files
    raw_bed_ref_child = convert_to_bedtool(raw_bed_child, sort=False)
    raw_bed_ref_parent = convert_to_bedtool(raw_bed_parent, sort=False)
    raw_bed_ref_depth_child = convert_to_bedtool(raw_bed_depth_child, sort=False)
    raw_bed_ref_depth_parent = convert_to_bedtool(raw_bed_ref_depth_parent_subset, sort=False)

    # Reformat de novo filt calls
    bed_child['chrom_type_sample'] = bed_child['chrom'] + "_" + bed_child['SVTYPE'] + "_" + bed_child['sample']
    bed_child['chrom_type_family'] = bed_child['chrom'] + "_" + bed_child['SVTYPE'] + "_" + bed_child['family_id'].astype(str)

    cols_keep_child = ['chrom_type_sample', 'start', 'end', 'name', 'svtype', 'sample', 'name_famid']
    cols_keep_parent = ['chrom_type_family', 'start', 'end', 'name', 'svtype', 'sample', 'name_famid']

    # Subset bed file down to just those where is_de_novo is still = True to make searching faster
    bed_child_de_novo = bed_child[bed_child['is_de_novo']]

    # INSERTIONS:
    verbose_print('Checking insertions in raw files', verbose)
    start = time.time()
    bed_filt_ins = bed_child_de_novo[bed_child_de_novo['SVTYPE'] == 'INS']
    if len(bed_filt_ins.index) > 0:
        verbose_print('Checking if insertion in proband is in raw files', verbose)
        bed_filt_ins_proband = convert_to_bedtool(bed_filt_ins, cols_to_keep=cols_keep_child, sort=True)
        ins_names_overlap_proband = get_insertion_intersection(bed_filt_ins_proband, raw_bed_ref_child)
        verbose_print('Checking if insertion in proband are also in raw files for the parents', verbose)
        bed_filt_ins_fam = convert_to_bedtool(bed_filt_ins, cols_to_keep=cols_keep_parent, sort=True)
        ins_names_overlap_parent = get_insertion_intersection(bed_filt_ins_fam, raw_bed_ref_parent)
        ins_names_overlap = [x for x in ins_names_overlap_proband if x not in ins_names_overlap_parent]
    else:
        ins_names_overlap = ['']
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Large CNVS: Reciprocal overlap large_raw_overlap
    verbose_print('Checking large cnvs in raw files', verbose)
    start = time.time()

    large_bed_filt_cnv_depth = bed_child_de_novo[(bed_child_de_novo['is_large_cnv']) & (bed_child_de_novo['SVLEN'] >= 5000)]
    large_bed_filt_cnv_other = bed_child_de_novo[(bed_child_de_novo['is_large_cnv']) & (bed_child_de_novo['SVLEN'] < 5000)]

    if len(large_bed_filt_cnv_other.index) > 0:
        verbose_print('Checking if intermediate cnv in proband is in raw files', verbose)
        bed_filt_cnv_proband_other = convert_to_bedtool(large_bed_filt_cnv_other, cols_to_keep=cols_keep_child, sort=True)
        large_cnv_names_overlap_proband_other = get_cnv_intersection_other(bed_filt_cnv_proband_other, raw_bed_ref_child, large_raw_overlap)
        verbose_print('Checking if intermediate cnvs in proband are also in raw files for the parents', verbose)
        bed_filt_cnv_fam_other = convert_to_bedtool(large_bed_filt_cnv_other, cols_to_keep=cols_keep_parent, sort=True)
        large_cnv_names_overlap_parent_other = get_cnv_intersection_other(bed_filt_cnv_fam_other, raw_bed_ref_parent, large_raw_overlap)
        large_cnv_names_overlap_other = [x for x in large_cnv_names_overlap_proband_other if x not in large_cnv_names_overlap_parent_other]
    else:
        large_cnv_names_overlap_other = ['']

    if len(large_bed_filt_cnv_depth.index) > 0:
        verbose_print('Checking if large cnv in proband is in raw files', verbose)
        bed_filt_cnv_proband_depth = convert_to_bedtool(large_bed_filt_cnv_depth, cols_to_keep=cols_keep_child, sort=True)
        large_cnv_names_overlap_proband_depth = get_cnv_intersection_depth(bed_filt_cnv_proband_depth, raw_bed_ref_depth_child, large_raw_overlap)
        verbose_print('Checking if large cnv in proband are also in raw files for the parents', verbose)
        bed_filt_cnv_fam_depth = convert_to_bedtool(large_bed_filt_cnv_depth, cols_to_keep=cols_keep_parent, sort=True)
        large_cnv_names_overlap_parent_depth = get_cnv_intersection_depth(bed_filt_cnv_fam_depth, raw_bed_ref_depth_parent, large_raw_overlap)
        large_cnv_names_overlap_depth = [x for x in large_cnv_names_overlap_proband_depth if x not in large_cnv_names_overlap_parent_depth]
    else:
        large_cnv_names_overlap_depth = ['']

    large_cnv_names_overlap = large_cnv_names_overlap_other + large_cnv_names_overlap_depth

    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Small CNVs - Reciprocal overlap small_raw_overlap
    verbose_print('Checking small cnvs in raw files', verbose)
    start = time.time()
    small_bed_filt_cnv = bed_child_de_novo[(bed_child_de_novo['is_small_cnv'])]  # We do not want to check against raw depth files if CNV <= 5kb
    if len(small_bed_filt_cnv.index) > 0:
        verbose_print('Checking if small cnv in proband is in raw files', verbose)
        bed_filt_cnv_proband = convert_to_bedtool(small_bed_filt_cnv, cols_to_keep=cols_keep_child, sort=True)
        small_cnv_names_overlap_proband = get_cnv_intersection_other(bed_filt_cnv_proband, raw_bed_ref_child, small_raw_overlap)
        verbose_print('Checking small cnvs in probands are also in raw files for parents', verbose)
        bed_filt_cnv_fam = convert_to_bedtool(small_bed_filt_cnv, cols_to_keep=cols_keep_parent, sort=True)
        small_cnv_names_overlap_parent = get_cnv_intersection_other(bed_filt_cnv_fam, raw_bed_ref_parent, small_raw_overlap)
        small_cnv_names_overlap = [x for x in small_cnv_names_overlap_proband if x not in small_cnv_names_overlap_parent]
    else:
        small_cnv_names_overlap = ['']
    delta = end - start
    print("Took %f seconds to process" % delta)

    bed_child.loc[~(bed_child['name_famid'].isin(ins_names_overlap + large_cnv_names_overlap + small_cnv_names_overlap)) & bed_child['is_de_novo'], 'filter_flag'] = 'not_in_raw_files_or_in_parents_raw_files'
    bed_child.loc[~(bed_child['name_famid'].isin(ins_names_overlap + large_cnv_names_overlap + small_cnv_names_overlap)) & bed_child['is_de_novo'], 'is_de_novo'] = False

    #############
    # FILTERING #
    #############
    verbose_print('Filtering out calls', verbose)

    # 1. Filter out calls in exclude regions
    verbose_print('Filtering out calls in exclude regions', verbose)
    start = time.time()
    # Reformat exclude_regions to bedtool
    exclue_regions_bt = convert_to_bedtool(exclude_regions, sort=True)
    # convert bed_final to bedtool
    cols_keep_exclude_regions = ['chrom', 'start', 'end', 'name', 'svtype', 'sample', 'name_famid']
    if len(bed_child.index) > 0:
        bed_child_bt = convert_to_bedtool(bed_child, cols_keep_exclude_regions, sort=True)
        exclude_regions_intersect = bed_child_bt.coverage(exclue_regions_bt).to_dataframe(disable_auto_names=True, header=None)  # HB said to use bedtools coverage, -f and -F give the same SVs to be removed
        if len(exclude_regions_intersect) != 0:
            remove_regions = exclude_regions_intersect[exclude_regions_intersect[10] > 0.5][6].to_list()
        else:
            remove_regions = ['']
    else:
        remove_regions = ['']

    bed_child.loc[bed_child['name_famid'].isin(remove_regions) & bed_child['is_de_novo'], 'filter_flag'] = 'in_blacklist_region'
    bed_child.loc[bed_child['name_famid'].isin(remove_regions) & bed_child['is_de_novo'], 'is_de_novo'] = False

    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # 2. Filter by type
    # Keep any SV type that is not a DEL, DUP, or INS
    keep_other_sv = bed_child[((~bed_child['SVTYPE'].isin(['DEL', 'DUP', 'INS'])) & (bed_child['filter_flag'] != 'AF') & (bed_child['filter_flag'] != 'in_parent'))]['name_famid'].to_list()

    # 3. Filter on DELs, DUPs, and INS
    # Filter by size
    # Filter out if large CNVs have parents overlap
    verbose_print('Filtering out large CNVs with parents overlap', verbose)
    start = time.time()
    remove_large = bed_child[(bed_child['is_large_cnv']) &
                             (bed_child['overlap_parent'])]['name_famid'].to_list()
    bed_child.loc[bed_child['name_famid'].isin(remove_large) & bed_child['is_de_novo'], 'filter_flag'] = 'large_cnv_with_parent_overlap'
    bed_child.loc[bed_child['name_famid'].isin(remove_large) & bed_child['is_de_novo'], 'is_de_novo'] = False
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    bed_child.to_csv(path_or_buf="after_large_cnv.txt", mode='a', index=False, sep='\t', header=True)

    # Filter out if small cnvs that are SR-only don't have BOTHSIDES_SUPPORT
    verbose_print('Filtering out small CNVs that are SR-only and dont have BOTHSIDES_SUPPORT', verbose)
    start = time.time()
    remove_small = bed_child[(bed_child['is_small_cnv']) &
                             (bed_child['EVIDENCE_FIX'] == 'SR') &
                             ~(bed_child.FILTER.str.contains('BOTHSIDES_SUPPORT'))]['name_famid'].to_list()
    bed_child.loc[bed_child['name_famid'].isin(remove_small) & bed_child['is_de_novo'], 'filter_flag'] = 'small_cnv_SR_only'
    bed_child.loc[bed_child['name_famid'].isin(remove_small) & bed_child['is_de_novo'], 'is_de_novo'] = False
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Filter out calls that are depth only and < depth_only_size
    verbose_print('Filtering out calls that are depth only and < depth_only_size', verbose)
    start = time.time()
    remove_depth_small = bed_child[(bed_child['is_depth_only_small'])]['name_famid'].to_list()
    bed_child.loc[bed_child['name_famid'].isin(remove_depth_small) & bed_child['is_de_novo'], 'filter_flag'] = 'depth_only_small_variant'
    bed_child.loc[bed_child['name_famid'].isin(remove_depth_small) & bed_child['is_de_novo'], 'is_de_novo'] = False
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Filter out DELs that are >500bp and RD_CN=2 and PE only envidence
    verbose_print('Filtering out DELs that are RD_CN=2 and PE only evidence', verbose)
    start = time.time()
    remove_dels = bed_child[(bed_child['SVTYPE'] == 'DEL') & ((bed_child['RD_CN'] == '2') | (bed_child['RD_CN'] == '3')) & (bed_child['EVIDENCE'] == 'PE')]['name_famid'].to_list()
    bed_child.loc[bed_child['name_famid'].isin(remove_dels) & bed_child['is_de_novo'], 'filter_flag'] = 'del_with_rdcn2_pe_only'
    bed_child.loc[bed_child['name_famid'].isin(remove_dels) & bed_child['is_de_novo'], 'is_de_novo'] = False
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # 4. Filter by quality
    # Filter out if parents GQ is <= gq_min
    verbose_print('Filtering if parents GQ <= min_gq', verbose)
    start = time.time()
    if len(bed_child.index) > 0:
        bed_child['keep_gq'] = bed_child.apply(lambda r: min_gq(r, gq_min), axis=1)
        remove_gq = bed_child[bed_child['keep_gq'] == 'Remove']['name_famid'].to_list()
        bed_child.loc[bed_child['name_famid'].isin(remove_gq) & bed_child['is_de_novo'], 'filter_flag'] = 'bad_gq_in_parents'
        bed_child.loc[bed_child['name_famid'].isin(remove_gq) & bed_child['is_de_novo'], 'is_de_novo'] = False
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Filter out INS that are dragen, manta or melt only and are SR only, have GQ=0, and FILTER contains 'HIGH_SR_BACKGROUND'
    verbose_print('Filtering out INS that are dragen, manta or melt only and SR only, with GQ=0 and FILTER contains HIGH_SR_BACKGROUND', verbose)
    start = time.time()
    remove_ins = bed_child[(bed_child['SVTYPE'] == 'INS') & ((bed_child['ALGORITHMS'] == 'dragen') | (bed_child['ALGORITHMS'] == 'manta') | (bed_child['ALGORITHMS'] == 'melt')) & (bed_child['EVIDENCE_FIX'] == 'SR') & ((bed_child['GQ'] == '0') | (bed_child.FILTER.str.contains('HIGH_SR_BACKGROUND')))]['name_famid'].to_list()
    bed_child.loc[bed_child['name_famid'].isin(remove_ins) & bed_child['is_de_novo'], 'filter_flag'] = 'ins_filter'
    bed_child.loc[bed_child['name_famid'].isin(remove_ins) & bed_child['is_de_novo'], 'is_de_novo'] = False
    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # Remove if low coverage in parents
    start = time.time()
    bed_child_coverage = bed_child[bed_child['is_de_novo']]
    verbose_print('Removing calls if there is low coverage evidence in parents', verbose)
    if len(bed_child_coverage.index) > 0:
        bed_child_coverage['median_coverage'] = bed_child_coverage.apply(
            lambda r: get_median_coverage(
                find_coverage(r, ped, sample_batches, bincov, family_member='mother'),
                find_coverage(r, ped, sample_batches, bincov, family_member='father'),
                coverage_cutoff),
            axis=1)
        remove_coverage = bed_child_coverage[bed_child_coverage['median_coverage'] == 'Remove']['name_famid'].to_list()
    else:
        remove_coverage = ['']

    bed_child.loc[bed_child['name_famid'].isin(remove_coverage) & bed_child['is_de_novo'], 'filter_flag'] = 'low_coverage_in_parents'
    bed_child.loc[bed_child['name_famid'].isin(remove_coverage) & bed_child['is_de_novo'], 'is_de_novo'] = False

    end = time.time()
    delta = end - start
    print("Took %f seconds to process" % delta)

    # 5. Clean up and remove duplicated CPX SV
    # Keep SVs
    bed_child.loc[bed_child['name_famid'].isin(keep_other_sv), 'is_de_novo'] = True
    bed_child.loc[bed_child['name_famid'].isin(keep_other_sv), 'filter_flag'] = 'not_del_dup_ins'

    # Remove duplicated CPX events that come from vcf output of module 18
    bed_child['is_cpx'] = (bed_child['SVTYPE'] == "CPX")
    bed_child['is_duplicated'] = bed_child.duplicated(keep='first', subset=['start', 'end', 'sample'])
    remove_duplicated = bed_child[(bed_child['is_duplicated']) & (bed_child['is_cpx'])]['name_famid'].to_list()
    bed_child.loc[bed_child['name_famid'].isin(remove_duplicated) & bed_child['is_de_novo'], 'is_de_novo'] = False
    bed_child.loc[bed_child['name_famid'].isin(remove_duplicated) & bed_child['is_de_novo'], 'filter_flag'] = 'duplicated_cpx'

    # Move sample column to column 6
    bed_final = bed_child
    col = bed_final.pop("sample")
    bed_final.insert(loc=5, column='sample', value=col, allow_duplicates=True)

    # Define output files
    output = bed_final
    de_novo = bed_final[(bed_final['is_de_novo']) | (bed_final['filter_flag'] == 'ins_filter') | (bed_final['in_gd'])]

    # Write output
    output.to_csv(path_or_buf=out_file, mode='a', index=False, sep='\t', header=True)
    de_novo.to_csv(path_or_buf=de_novo_out_file, mode='a', index=False, sep='\t', header=True)


if __name__ == '__main__':
    main()
