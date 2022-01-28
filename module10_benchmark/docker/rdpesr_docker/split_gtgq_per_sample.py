import os
import gzip
import pysam
import argparse

def sample_list_readin(sample_name):
	sample_list = []
	fin=open(sample_name)
	for line in fin:
		pin=line.strip().split()
		if pin[0][:2]=="__":
			sample_list.append(pin[0])
	fin.close()
	return sample_list

def split_gtgq_per_sample(sample_list, vcf_file, chr_name):
	fin=pysam.VariantFile(vcf_file)
	fout = {sample: gzip.open("per_sample_GTGQ/"+sample+'.'+chr_name+'.gtgq.gz', mode='w', compresslevel=3) for sample in sample_list}
	for sample in sample_list:
		out_str = '\t'.join(['SVID','GT', 'CN', 'CNQ', 'EV', 'GQ', 'PE_GQ', 'PE_GT', 'RD_CN', 'RD_GQ', 'SR_GQ', 'SR_GT']) + '\n'
		fout[sample].write(out_str.encode())

	for record in fin:
		print(record.id)
		for sample in sample_list:
			sample_data = record.samples[sample]
			if sample_data['GT']==(0,0) or sample_data['GT']==(None,None): 
				continue
			out_str = '\t'.join([str(j) for j in [record.id, '/'.join([str(i) for i in sample_data['GT']]),sample_data['CN'],sample_data['CNQ'],
										','.join([i for i in sample_data['EV']]),sample_data['GQ'],
										sample_data['RD_CN'],sample_data['RD_GQ'],
										sample_data['PE_GT'],sample_data['PE_GQ'],
										sample_data['SR_GT'],sample_data['SR_GQ']]]) + '\n'
			fout[sample].write(out_str.encode())

def main():
	parser = argparse.ArgumentParser(description='Split gtgq info per sample')
	parser.add_argument('sample', metavar='', type=str,
						help='a 1-column list with names of samples to process')
	parser.add_argument('vcf', metavar='', type=str,
						help='vcf to parse from')
	parser.add_argument('chr', metavar='', type=str,
						help='name of chromoosome to process')

	args = parser.parse_args()
	sample_list = sample_list_readin(args.sample)
	split_gtgq_per_sample(sample_list, args.vcf, args.chr)

if __name__ == "__main__":
	main()



