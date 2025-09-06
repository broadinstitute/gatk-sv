#!/usr/bin/env python3

import argparse
import csv
import gzip
import sys
import re
from pathlib import Path

def parse_info_field(info_field):
    """Parse VCF INFO field into dictionary"""
    info_dict = {}
    if info_field == ".":
        return info_dict
    
    for item in info_field.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    return info_dict

def parse_filter_field(filter_field):
    """Parse VCF FILTER field into list"""
    if filter_field == "." or filter_field == "PASS":
        return [filter_field]
    return filter_field.split(",")

def load_reference_genes(ref_genes_file):
    """Load reference genes file into dictionary"""
    genes_dict = {}
    if not ref_genes_file:
        return genes_dict
        
    with open(ref_genes_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene_name = row['gene_name']
            genes_dict[gene_name] = {
                'LOEUF': float(row['LOEUF']) if row['LOEUF'] != '' else None,
                'LOEUF_tile': int(row['LOEUF_tile']) if row['LOEUF_tile'] != '' else None
            }
    return genes_dict

def check_gene_criteria(info_dict, infos_list, genes_dict, min_loeuf=None, max_loeuf=None, min_loeuf_tile=None, max_loeuf_tile=None):
    """Check if variant meets gene-based criteria"""
    if not infos_list or not genes_dict:
        return True
        
    for info_key in infos_list:
        if info_key in info_dict:
            gene_names = info_dict[info_key].split(",")
            for gene_name in gene_names:
                if gene_name in genes_dict:
                    gene_data = genes_dict[gene_name]
                    
                    # Check LOEUF criteria
                    if min_loeuf is not None or max_loeuf is not None:
                        loeuf_val = gene_data.get('LOEUF')
                        if loeuf_val is not None:
                            if min_loeuf is not None and loeuf_val < min_loeuf:
                                continue
                            if max_loeuf is not None and loeuf_val >= max_loeuf:
                                continue
                            return True
                    
                    # Check LOEUF_tile criteria
                    if min_loeuf_tile is not None or max_loeuf_tile is not None:
                        loeuf_tile_val = gene_data.get('LOEUF_tile')
                        if loeuf_tile_val is not None:
                            if min_loeuf_tile is not None and loeuf_tile_val < min_loeuf_tile:
                                continue
                            if max_loeuf_tile is not None and loeuf_tile_val >= max_loeuf_tile:
                                continue
                            return True
    
    return False

def meets_criteria(variant_fields, stratification, genes_dict):
    """Check if variant meets all criteria for a stratification"""
    chrom, pos, var_id, ref, alt, qual, filter_field, info_field = variant_fields[:8]
    
    info_dict = parse_info_field(info_field)
    filter_list = parse_filter_field(filter_field)
    
    # Check SV types
    if stratification.get('sv_types'):
        svtype = info_dict.get('SVTYPE')
        if not svtype or svtype not in stratification['sv_types']:
            return False
    
    # Check size criteria
    svlen = info_dict.get('SVLEN')
    if svlen:
        try:
            svlen_val = abs(int(svlen))
            if stratification.get('min_size') is not None and svlen_val < stratification['min_size']:
                return False
            if stratification.get('max_size') is not None and svlen_val >= stratification['max_size']:
                return False
        except ValueError:
            pass
    
    # Check filter criteria
    if stratification.get('filters'):
        if not any(f in stratification['filters'] for f in filter_list):
            return False
        if any(f not in stratification['filters'] for f in filter_list if f != '.'):
            return False
    
    # Check info criteria
    if stratification.get('infos'):
        if not any(info_key in info_dict for info_key in stratification['infos']):
            return False
    
    # Check AC/AF criteria
    for field, min_key, max_key in [('AC', 'min_ac', 'max_ac'), ('AF', 'min_af', 'max_af')]:
        if field in info_dict:
            try:
                val = float(info_dict[field])
                if stratification.get(min_key) is not None and val < stratification[min_key]:
                    return False
                if stratification.get(max_key) is not None and val >= stratification[max_key]:
                    return False
            except ValueError:
                pass
    
    # Check gene criteria
    if not check_gene_criteria(
        info_dict, 
        stratification.get('infos'),
        genes_dict,
        stratification.get('min_loeuf'),
        stratification.get('max_loeuf'),
        stratification.get('min_loeuf_tile'),
        stratification.get('max_loeuf_tile')
    ):
        return False
    
    return True

def parse_stratifications(strat_file):
    """Parse plot stratifications file"""
    stratifications = []
    plot_names = set()
    
    with open(strat_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Validate plot_name
            plot_name = row['plot_name']
            if ' ' in plot_name or plot_name in plot_names:
                raise ValueError(f"Invalid plot_name: {plot_name}")
            plot_names.add(plot_name)
            
            # Parse stratification
            strat = {'plot_name': plot_name}
            
            # Parse optional fields
            for field, default_val in [
                ('sv_types', None), ('filters', None), ('infos', None),
                ('min_size', None), ('max_size', None),
                ('min_ac', None), ('max_ac', None),
                ('min_af', None), ('max_af', None),
                ('min_loeuf', None), ('max_loeuf', None),
                ('min_loeuf_tile', None), ('max_loeuf_tile', None)
            ]:
                val = row.get(field, '').strip()
                if val and val != '-1':
                    if field in ['sv_types', 'filters', 'infos']:
                        strat[field] = [x.strip() for x in val.split(',')]
                    elif field.startswith('min_') or field.startswith('max_'):
                        try:
                            strat[field] = float(val) if 'af' in field or 'loeuf' in field else int(val)
                        except ValueError:
                            pass
            
            # Validate plot type
            depth_plot = row.get('depth_plot', '').strip().lower() == 'true'
            igv_plot = row.get('igv_plot', '').strip().lower() == 'true'
            
            if depth_plot and igv_plot:
                raise ValueError(f"Both depth_plot and igv_plot are True for {plot_name}")
            if not depth_plot and not igv_plot:
                raise ValueError(f"Neither depth_plot nor igv_plot is True for {plot_name}")
            
            strat['depth_plot'] = depth_plot
            strat['igv_plot'] = igv_plot
            
            stratifications.append(strat)
    
    return stratifications

def main():
    parser = argparse.ArgumentParser(description='Filter VCF based on plot stratifications')
    parser.add_argument('vcf_file', help='Input VCF file (can be gzipped)')
    parser.add_argument('stratifications_file', help='Plot stratifications TSV file')
    parser.add_argument('output_prefix', help='Output prefix for filtered VCF files')
    parser.add_argument('--reference-genes', help='Reference genes TSV file')
    
    args = parser.parse_args()
    
    # Load reference genes if provided
    genes_dict = load_reference_genes(args.reference_genes) if args.reference_genes else {}
    
    # Parse stratifications
    stratifications = parse_stratifications(args.stratifications_file)
    
    # Open input VCF
    if args.vcf_file.endswith('.gz'):
        vcf_file = gzip.open(args.vcf_file, 'rt')
    else:
        vcf_file = open(args.vcf_file, 'r')
    
    # Initialize output files
    output_files = {}
    for strat in stratifications:
        plot_name = strat['plot_name']
        output_file = f"{args.output_prefix}.{plot_name}.vcf"
        output_files[plot_name] = open(output_file, 'w')
    
    # Process VCF
    header_written = {plot_name: False for plot_name in output_files.keys()}
    
    for line in vcf_file:
        if line.startswith('#'):
            # Write header to all output files
            for f in output_files.values():
                f.write(line)
            for plot_name in header_written:
                header_written[plot_name] = True
        else:
            # Process variant
            fields = line.strip().split('\t')
            for strat in stratifications:
                if meets_criteria(fields, strat, genes_dict):
                    output_files[strat['plot_name']].write(line)
    
    # Close files
    vcf_file.close()
    for f in output_files.values():
        f.close()
    
    # Create summary of output files
    summary_data = []
    for strat in stratifications:
        plot_name = strat['plot_name']
        output_file = f"{args.output_prefix}.{plot_name}.vcf"
        summary_data.append({
            'plot_name': plot_name,
            'depth_plot': strat['depth_plot'],
            'igv_plot': strat['igv_plot'],
            'vcf_file': output_file
        })
    
    # Write summary file
    summary_file = f"{args.output_prefix}.summary.tsv"
    with open(summary_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=['plot_name', 'depth_plot', 'igv_plot', 'vcf_file'], delimiter='\t')
        writer.writeheader()
        writer.writerows(summary_data)
    
    print(f"Created {len(stratifications)} filtered VCF files")
    print(f"Summary written to {summary_file}")

if __name__ == '__main__':
    main() 