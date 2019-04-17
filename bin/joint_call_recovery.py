import sys
import re

annotated_merged_f = open(sys.argv[1])
LD_dels_support_reads_threshold = int(sys.argv[2])
SID_dels_number_threshold = int(sys.argv[3])
LD_dels_number_threshold = int(sys.argv[4])
LD_sr_number_threshold = int(sys.argv[5])
del_len_threshold = int(sys.argv[6])

def LD_dels_LD_sr_dels_eligibility(LD_dels, LD_sr_dels):
	if LD_dels >= LD_dels_number_threshold and LD_sr_dels >= LD_sr_number_threshold:
		dele_stat = "Rescue"
	else:
		dele_stat = "Exclude"
	return dele_stat
	
for line in annotated_merged_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "chr":
		print("Deletion_status", "\t".join(line_l[0:11]), "Samples_num", "Deletion_classes", "\t".join(line_l[11:]), sep="\t")
		continue
	dele_stat = "Exclude_B"
	acen_line = line_l[6]
	samples_line = line_l[13:]
	samples_list = []
	deletion_len = int(line_l[2]) - int(line_l[1])
	if deletion_len == 0:
		deletion_len = 1
	for sample in samples_line:
		if not sample == "-":
			samples_list.append(sample)
	num_samples = len(samples_list)
	SID_dels = 0
	F_dels = 0
	B_dels = 0
	LD_dels = 0
	LD_sr_dels = 0
	LI_dels = 0
	for sample in samples_list:
		sample_l = re.split(",", sample)
		for samp in sample_l:
			samp_l = re.split(";", samp)
			if len(samp_l) == 1:
				continue
			del_info = samp_l[9]
			del_l = re.split("_", del_info)
			del_class = del_l[1]
			del_supp_reads = int(del_l[2])
			if del_class == "SID":
				SID_dels += 1
			if del_class == "B":
				B_dels += 1
			if del_class == "F":
				F_dels += 1
			if del_class == "LD":
				LD_dels += 1
				if del_supp_reads >= LD_dels_support_reads_threshold:
					LD_sr_dels += 1
			if del_class == "LI":
				LI_dels += 1
	dele_class = "SID:" + str(SID_dels) + ";" + "B:" + str(B_dels) + ";" + "F:" + str(F_dels) + ";" + "LD:" + str(LD_dels) + ";" + "LD_sr:" + str(LD_sr_dels) + ";" + "LI:" + str(LI_dels)
	res_line = ["\t".join(line_l[0:13]), str(num_samples), dele_class, "\t".join(line_l[13:])]
	res_line_str = "\t".join(res_line)
	if deletion_len < del_len_threshold:
		if LI_dels > 0:
			dele_stat = "Exclude_LI"
		elif SID_dels != 0 and SID_dels >= SID_dels_number_threshold:
			dele_stat = "Rescue"
		elif SID_dels != 0 and SID_dels < SID_dels_number_threshold:
			dele_stat = "Exclude"
		elif B_dels != 0 and F_dels != 0:
			dele_stat = "Rescue"
		elif (SID_dels != 0 and B_dels != 0) or (SID_dels != 0 and F_dels != 0):
			dele_stat = "Rescue"
		elif (B_dels != 0 or F_dels != 0 or SID_dels != 0) and LD_dels != 0:
			dele_stat = "Rescue"
		elif B_dels == 0 and F_dels == 0 and SID_dels == 0 and LI_dels == 0 and LD_dels != 0:
			dele_stat = LD_dels_LD_sr_dels_eligibility(LD_dels, LD_sr_dels)
	elif deletion_len >= del_len_threshold:
		if LI_dels >0:
			dele_stat = "Exclude_LI"
		elif SID_dels != 0 and B_dels == 0 and F_dels == 0 and LD_dels == 0 and LD_sr_dels == 0:
			dele_stat = "Exclude"
		elif B_dels != 0 and F_dels != 0:
			dele_stat = "Rescue" 
		elif (B_dels != 0 or F_dels != 0) and (SID_dels != 0 and LD_dels == 0):
			dele_stat = "Exclude_a"
		elif (B_dels != 0 or F_dels != 0 or SID_dels != 0) and LD_dels != 0:
			dele_stat = LD_dels_LD_sr_dels_eligibility(LD_dels, LD_sr_dels)
		elif B_dels != 0 and F_dels != 0 and SID_dels != 0 and LD_dels != 0:
			dele_stat = LD_dels_LD_sr_dels_eligibility(LD_dels, LD_sr_dels)
		elif B_dels == 0 and F_dels == 0 and SID_dels == 0 and LI_dels == 0 and LD_dels != 0:
				dele_stat = LD_dels_LD_sr_dels_eligibility(LD_dels, LD_sr_dels)
	if not acen_line == "-":
		dele_stat = "Acen"
	print(dele_stat, res_line_str, sep="\t")
