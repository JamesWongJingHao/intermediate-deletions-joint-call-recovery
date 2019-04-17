import sys
import re

HWE_f = open(sys.argv[1])
del_res_f = open(sys.argv[2])

dele_HWE_dict = {}
for line in HWE_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	line_l_s = []
	for entry in line_l:
		entry2 = entry.replace("\"", '')
		line_l_s.append(entry2)
	if line_l_s[0] == "obs_homo":
		continue
	deletion = line_l_s[0] 
	p_val = line_l[-1]
	if p_val == "1":
		p_val = "1.0"
	if not deletion in dele_HWE_dict:
		dele_HWE_dict[deletion] = p_val 
for row in del_res_f:
	row = row.replace('\n', '')
	row_l = re.split("\t", row)
	if row_l[0] == "chr":
		print("\t".join(row_l[0:14]), "obs_Homo", "obs_hetero", "obs_HomoWT", "exp_Homo", "exp_Hetero", "exp_HomoWT", "P_HWE", "\t".join(row_l[14:]), sep="\t")
		continue
	deletion_name = row_l[0] + "_" + row_l[1] + "_" + row_l[2]
#	samples = row_l[16:]
#	number_samples = len(samples)
	samples = []
	samples_zygosity = []
	sampleline = row_l[16:]
	number_samples = len(sampleline)
	for sample in sampleline:
		if not sample == "-":
			samples.append(sample)
	zygosity_list = []
	for sample in samples:
		sample_l = re.split(",", sample)
		if len(sample_l) == 1:
			sampl = re.split(";", sample)
			zygosity_list.append(sampl[2])
		elif len(sample_l) > 1:
			sampl = sample_l[0]
			samp = re.split(";", sample)
			samp_l = samp[2]
			zygosity_list.append(samp_l)
	obs_Hete_cnt = 0
	obs_Homo_cnt = 0
	Skip_cnt = 0
	for zygos in zygosity_list:
		if zygos == "Homo":
			obs_Homo_cnt += 1
		elif zygos == "Hete":
			obs_Hete_cnt += 1
		elif zygos == "SKIP":
			Skip_cnt + 1
	obs_WT_homo = number_samples - (obs_Homo_cnt + obs_Hete_cnt)
	total = obs_Homo_cnt + obs_Hete_cnt + obs_WT_homo
	P_del = ((2*obs_Homo_cnt) + obs_Hete_cnt)/(total*2)
	P_WT = ((2*obs_WT_homo) + obs_Hete_cnt)/(total*2)
	exp_homo = round(((P_del*P_del)*total), 0)
	exp_hete = round((((2*P_del)*(1 - P_del))*total), 0)
	exp_homo_WT = round((((1 - P_del)*(1 - P_del))*total), 0)
	if deletion_name in dele_HWE_dict:
		print("\t".join(row_l[0:16]), obs_Homo_cnt, obs_Hete_cnt, obs_WT_homo, str(int(exp_homo)), str(int(exp_hete)), str(int(exp_homo_WT)),  dele_HWE_dict[deletion_name], "\t".join(row_l[16:]), sep="\t") 
