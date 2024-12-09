#!python script to calculate Hudson Fst per SV site across all the populations
import allel
import pysam
import argparse
import numpy as np
import pandas as pd
from itertools import combinations

# Hudson’s Fst calculation
def hudson_fst(ac1, ac2):
    p1 = ac1[1] / ac1[0] if ac1[0] > 0 else 0
    p2 = ac2[1] / ac2[0] if ac2[0] > 0 else 0
    num = (p1 - p2) ** 2 - (p1 * (1 - p1) / (ac1[0] - 1) if ac1[0] > 1 else 0) - (p2 * (1 - p2) / (ac2[0] - 1) if ac2[0] >1 else 0)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return num, den

#filter for common variants
def filter_common_variants(record, pop_indices, min_af=0.01, min_samples=1000):
    # Ensure global AF > 0.01
    global_af = record.info['AF'][0]
    if global_af <= min_af:
        return False
    # Ensure AF > 0.01 in each population with >1000 samples
    for pop in pop_indices.keys():
        if len(pop_indices[pop]) >= min_samples:
            pop_af = record.info[pop+'_AF'][0]
            if pop_af <= min_af:
                return False
    return True

#function to calculate fst table for each pair of populations for each loci
def calculate_fst_table(vcf,sample_pop):
        # File paths
        vcf_file = vcf
        pop_file = sample_pop  # Replace with your population file path
        # Load the population data
        pop_data = pd.read_csv(pop_file, sep='\t', header=None)
        pop_data.columns = ['sample_id','population']
        pop_dict = {pop: pop_data.loc[pop_data['population'] == pop, 'sample_id'].values for pop in pop_data['population'].unique()}
        pop_dict.pop("rmi", None)  # 'None' is the default value if "rmi" is not found
        # Open VCF using pysam
        vcf = pysam.VariantFile(vcf_file)
        # Verify samples in VCF are in population file
        vcf_samples = list(vcf.header.samples)
        pop_samples = pop_data['sample_id'].tolist()
        assert set(vcf_samples).issubset(set(pop_samples)), "Mismatch between VCF samples and population file sample IDs"
        # Map sample indices in VCF to populations
        pop_indices = {pop: [vcf_samples.index(s) for s in samples if s in vcf_samples] for pop, samples in pop_dict.items()}
        # Create a DataFrame to store num and den for each population pair at each site
        pop_pairs = list(combinations(pop_indices.keys(), 2))
        columns = ['ID'] + [f'{pop1}_{pop2}_num' for pop1, pop2 in pop_pairs] + [f'{pop1}_{pop2}_den' for pop1, pop2 in pop_pairs]
        fst_data = pd.DataFrame(columns=columns)
        # Iterate over each variant in the VCF
        for record in vcf:
            if not 'PASS' in record.filter.keys(): continue
            # Initialize row data
            row = {'ID': record.chrom+'_'+str(record.pos)}
            if filter_common_variants(record, pop_indices,0.01):
                    print(record.chrom+'_'+str(record.pos))
            # Calculate num and den for each population pair
                    for (pop1, pop2) in pop_pairs:
                        # Calculate allele counts for each population
                        ac1 = [record.info['AN_'+pop1], record.info['AC_'+pop1][0]]
                        ac2 = [record.info['AN_'+pop2], record.info['AC_'+pop2][0]]
                        # Calculate Hudson's Fst numerator and denominator
                        num, den = hudson_fst(ac1, ac2)
                        # Store in the row with appropriate column names
                        row[f'{pop1}_{pop2}_num'] = num
                        row[f'{pop1}_{pop2}_den'] = den
                    # Append row data to DataFrame
                    fst_data = fst_data._append(row, ignore_index=True)
        #integrate Fst across populations for each variant site:
        num_cols = [col for col in fst_data.columns if 'num' in col]
        den_cols = [col for col in fst_data.columns if 'den' in col]
        fst_data['sum_num'] = fst_data[num_cols].sum(axis=1)
        fst_data['sum_den'] = fst_data[den_cols].sum(axis=1)
        fst_data['fst'] = fst_data['sum_num']/fst_data['sum_den']
        return fst_data

#function to integrate fst across all SV sites
def integrate_fst_across_sites(fst_data):
        # Identify all unique population pairs from the column names
        population_pairs = set(col.replace('_num','') for col in fst_data.columns if 'num' in col)
        # Initialize a dictionary to store aggregated numerator and denominator for each population pair
        aggregated_fst = {pair: {'total_num': 0, 'total_den': 0} for pair in population_pairs}
        # Aggregate numerator and denominator across sites for each population pair
        for pair in population_pairs:
            num_col = f'{pair}_num'
            den_col = f'{pair}_den'
            # Sum numerators and denominators across all sites for the given population pair
            aggregated_fst[pair]['total_num'] = fst_data[num_col].sum()
            aggregated_fst[pair]['total_den'] = fst_data[den_col].sum()
        # Calculate Fst as total_num / total_den for each population pair and save results
        fst_aggregated_results = []
        for pair, values in aggregated_fst.items():
            if values['total_den'] > 0:
                fst_value = values['total_num'] / values['total_den']
            else:
                fst_value = None  # Set to None if denominator is zero
            fst_aggregated_results.append({
                'Population_Pair': pair,
                'Aggregated_Fst': fst_value,
                'Total_Num': values['total_num'],
                'Total_Den': values['total_den']
            })
        # Convert to DataFrame and save to CSV
        fst_aggregated_df = pd.DataFrame(fst_aggregated_results)
        return fst_aggregated_df

#define main function
def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Process a CSV file and calculate ratios")
    parser.add_argument("-v", "--vcf", required=False, help="input vcf that include variant sites and allele frequency annotations")
    parser.add_argument("-s", "--sample_pop", required=False, help="input file including sample ID and populations deliminated by tab")
    parser.add_argument("-o", "--output", required=False, help="output file that include per-site Fst information")
    parser.add_argument("-p", "--population_output", required=False, help="output file that include Fst for each pair of populations")

    # Parse arguments
    args = parser.parse_args()

    fst_data = calculate_fst_table(args.vcf, args.sample_pop)

    # Save per-site Fst components (num and den) to CSV
    fst_data.to_csv(args.output, index=False, sep='\t')

    #calculate and save Fst of per pair of populations
    fst_aggregated_df = integrate_fst_across_sites(fst_data)
    fst_aggregated_df.to_csv(args.population_output, index=False)


if __name__ == "__main__":
    main()
