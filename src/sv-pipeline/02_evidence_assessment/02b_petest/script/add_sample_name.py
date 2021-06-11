import argparse
parser = argparse.ArgumentParser("add_sample_name")
parser.add_argument("sample", type=str,
                    help="namd of sample to be extracted from vcf and kept in bed format")
parser.add_argument("-v", '--version', type=str,
                    help="indicate the version (old / new) of the vcf set")

# define workding and file locations:
args = parser.parse_args()
sample_name = args.sample
fin = open(sample_name)
fo = open(sample_name + '.with_name', 'w')
for line in fin:
    pin = line.strip().split()
    pin_new = pin + [sample_name.split('.')[0]]
    print('\t'.join(pin_new), file=fo)

fin.close()
fo.close()
