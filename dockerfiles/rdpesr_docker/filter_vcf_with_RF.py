import os
import sys
import pysam
import argparse

def readin_RF_cutoffs(cff_file):
  fin=open(cff_file)
  out = {}
  header = []
  for line in fin:
    pin=line.strip().split()
    if pin[0]=="svtype":
      header = pin
    else: 
      svtype = pin[header.index('svtype')]
      size_name = pin[header.index('size')]
      size_min = pin[header.index('size_min')]
      size_max = pin[header.index('size_max')]
      cff_pbsv = pin[header.index('cff_pbsv')]
      cff_vapor = pin[header.index('cff_vapor')]
      if not svtype in out.keys():
        out[svtype] = {}
      out[svtype][size_name] = [cff_pbsv, cff_vapor]
  return out

def decide_size_range(size):
  if not size>100:
    size_range = 'under100bp'
  else:
    if not size >1000:
      size_range = '100to1Kb'
    else:
      if not size>5000:
        size_range = '1to5Kb'
      else:
        size_range = 'over5Kb'
  return size_range

def main():
  parser = argparse.ArgumentParser(
      description=__doc__,
      formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('vcf_in', help='VCF to annotate.')
  parser.add_argument('vcf_out', help='Path to output VCF. Accepts "-" ' +
                      'and "stdout". Default: stdout', default='stdout')
  parser.add_argument('--stats', help='Optional path to output .tsv file ' +
                      'with NCR stats for all records in vcf_in.')
  args = parser.parse_args()

  cff_table = readin_RF_cutoffs(args.stats)

  fin=pysam.VariantFile(args.vcf_in)
  # Open connection to output VCF
  if args.vcf_out in '- stdout'.split():
      fo = pysam.VariantFile(sys.stdout, 'w', header=fin.header)
  else:
      fo = pysam.VariantFile(args.vcf_out, 'w', header=fin.header)

  for record in fin:
    print(record.id)
    if record.info['SVTYPE'] in cff_table.keys():
      if not '<INS:ME:ALU>' in record.alts and not '<INS:ME:LINE1>' in record.alts and not '<INS:ME:SVA>' in record.alts:
        size_range =  decide_size_range(record.info["SVLEN"])
        [cff_pbsv, cff_vapor] = cff_table[record.info['SVTYPE']][size_range]
      else:
        if '<INS:ME:ALU>' in record.alts:
          [cff_pbsv, cff_vapor] = cff_table['INS:ME:ALU']['all']
        if '<INS:ME:LINE1>' in record.alts:
          [cff_pbsv, cff_vapor] = cff_table['INS:ME:LINE1']['all']
        if '<INS:ME:SVA>' in record.alts:
          [cff_pbsv, cff_vapor] = cff_table['INS:ME:SVA']['all']
      cff_pbsv = float(cff_pbsv)
      cff_vapor = float(cff_vapor)
      for sample in record.samples.keys():
        if record.samples[sample]['GT']==(0,1) or record.samples[sample]['GT']==(1,1):
          if record.samples[sample]['PBRF']==None and record.samples[sample]['VaRF']==None:
            record.samples[sample]['GT']=(None, None)
          else:
            if record.samples[sample]['PBRF'] < cff_pbsv and record.samples[sample]['VaRF'] < cff_vapor:
              record.samples[sample]['GT']=(None, None)
      fo.write(record)

  fin.close()
  fo.close()

if __name__ == '__main__':
    main()


