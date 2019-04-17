import sys
import re
import pandas as pd

filtered_deletions_f = open(sys.argv[1])
number_samples = 0
total_obs_homo = 0
total_obs_hete = 0
for line in filtered_deletions_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "chr":
		samples = line_l[14:]
		number_samples = len(samples)
		print("Deletion", "obs_homo", "obs_hete", "obs_homoWT", "exp_homo", "exp_hete", "exp_homoWT", sep="\t")
		continue
	samples = []
	samples_zygosity = []
	sampleline = line_l[17:]	
	number_samples = len(sampleline)
	for sample in sampleline:
		if not sample == "-":
			samples.append(sample)
	deletion_name = line_l[0] + "_" + line_l[1] + "_" + line_l[2]
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
	obs_Homo_cnt = 0
	obs_Hete_cnt = 0
	Skip_cnt = 0
	for zygos in zygosity_list:
		if zygos == "Homo":
			obs_Homo_cnt += 1
		elif zygos == "Hete":
			obs_Hete_cnt += 1
		elif zygos == "SKIP":
			Skip_cnt + 1  
	total_obs_hete += obs_Hete_cnt
	obs_WT_homo = number_samples - (obs_Homo_cnt + obs_Hete_cnt)
	total = obs_Homo_cnt + obs_Hete_cnt + obs_WT_homo
	P_del = ((2*obs_Homo_cnt) + obs_Hete_cnt)/(total*2)
	P_WT = ((2*obs_WT_homo) + obs_Hete_cnt)/(total*2)
	exp_homo = round(((P_del*P_del)*total), 0)
	exp_hete = round((((2*P_del)*(1 - P_del))*total), 0)
	exp_homo_WT = round((((1 - P_del)*(1 - P_del))*total), 0)
	
	contingency = [deletion_name, str(obs_Homo_cnt), str(obs_Hete_cnt), str(obs_WT_homo), str(int(exp_homo)), str(int(exp_hete)), str(int(exp_homo_WT))]
	contingency_str = "\t".join(contingency)
	print(contingency_str)

