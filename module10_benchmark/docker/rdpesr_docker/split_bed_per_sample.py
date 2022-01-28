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

def split_bed_per_sample(sample_list, bed_file, chr_name):
	fin=os.popen(r'''zcat %s'''%(bed_file))
	fout = {sample: gzip.open("per_sample_bed/"+sample+'.'+chr_name+'.bed.gz', mode='w', compresslevel=3) for sample in sample_list}
	for sample in sample_list:
		out_str = '\t'.join(['#chrom','start','end','VID','svtype','length','AF','samples']) + '\n'
		fout[sample].write(out_str.encode())

	for line in fin:
		pin = line.strip().split()
		sample_data = pin[5].split(',')
		for sample in sample_list:
			if sample in sample_data:
				out_str = '\t'.join([str(j) for j in pin[:4]+pin[6:8]+[0,sample]]) + '\n'
				fout[sample].write(out_str.encode())

def main():
	parser = argparse.ArgumentParser(description='Split gtgq info per sample')
	parser.add_argument('sample', metavar='', type=str,
						help='a 1-column list with names of samples to process')
	parser.add_argument('bed', metavar='', type=str,
						help='bed to parse from')
	parser.add_argument('chr', metavar='', type=str,
						help='name of chromoosome to process')

	args = parser.parse_args()
	sample_list = sample_list_readin(args.sample)
	split_bed_per_sample(sample_list, args.bed, args.chr)

if __name__ == "__main__":
	main()




