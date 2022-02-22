#! python
import os
import argparse

def cff_table_readin(cff_table):
	fin=open(cff_table)
	out={}
	for line in fin:
		pin=line.strip().split()
		out[pin[0]]=float(pin[1])
	return out

def calculate_SVID_pass_rate(folder_path, boost_cutoff_table):
	SVID_pass_rate = {}
	for file in os.listdir(folder_path):
		if file.split('.')[1] == 'boost_filtered':
			file_name=folder_path + file
			fin=open(file_name)
			for line in fin:
				pin=line.strip().split()
				if pin[0]=='SVID': continue
				if not pin[0] in SVID_pass_rate.keys():
					SVID_pass_rate[pin[0]]=[0,0]
				svtype = pin[0].split('_')[2]
				boost_cutoff = boost_cutoff_table[svtype]
				if float(pin[1])>boost_cutoff:
					SVID_pass_rate[pin[0]][1]+=1
				else:
					SVID_pass_rate[pin[0]][0]+=1
			fin.close()
	return SVID_pass_rate

def write_output(out_name, SVID_pass_rate):
	fo=open(out_name,'w')
	print('\t'.join(['SVID','count_samp_fail','count_samp_pass']), file=fo)
	for i in SVID_pass_rate.keys():
		print('\t'.join([str(j) for j in [i]+SVID_pass_rate[i]]), file=fo)
	fo.close()

def main():
	parser = argparse.ArgumentParser(description='Integrate boost scores across samples')
	parser.add_argument('--cff_table', metavar='', type=str, help='table including svtype and corresponding boost score cutoff')
	parser.add_argument('--path', metavar='', type=str, help='name of path to input files')
	parser.add_argument('--output', metavar='', type=str, help='name of output file')
	args = parser.parse_args()
	boost_cutoff_table = cff_table_readin(args.cff_table)	
	folder_path = args.path
	if not folder_path[-1]=='/':
		folder_path+='/'
	out_name = args.output
	SVID_pass_rate = calculate_SVID_pass_rate(folder_path, boost_cutoff_table)
	write_output(out_name, SVID_pass_rate)

if __name__ == "__main__":
	main()


