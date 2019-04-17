import sys
import re

annotate_f = open(sys.argv[1])
acen_f = open(sys.argv[2])
simplerep_f = open(sys.argv[3])
repmask_f = open(sys.argv[4])
overall_res_f = open(sys.argv[5])
flank_range = int(sys.argv[6])
gene_dict = {}
for line in annotate_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "#bin":
		continue
	chrm = line_l[2]
	gene_dict.setdefault(chrm, []).append(line_l)
acen_dict = {}
for line in acen_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "chrom":
		continue
	chrom = line_l[0]
	acen = [line_l[1], line_l[2], line_l[6]]
	acen_dict.setdefault(chrom, []).append(acen)
simplerep_dict = {}
for line in simplerep_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "#bin":
		continue
	chrom = line_l[1]
	simrep = [line_l[2], line_l[3], line_l[4]]
	simplerep_dict.setdefault(chrom, []).append(simrep)
repmask_dict = {}
for line in repmask_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "bin":
		continue
	chrom = line_l[5]
	rep = [line_l[6], line_l[7], line_l[9], line_l[10], line_l[11], line_l[12]]
	repmask_dict.setdefault(chrom, []).append(rep)
output_str = ""
for row in overall_res_f:
	row = row.replace('\n', '')
	row_l = re.split("\t", row)
	if row_l[0] == "chr":
		print("\t".join(row_l[0:3]), "gene_name", "gene_annotation", "UTR", "centromere/telomere", "simple_repeats", "repeatmasker", "forward_flanking_region", "reverse_flanking_region", "\t".join(row_l[3:]), sep="\t")
		continue
	deletion_front_flank = int(row_l[1]) - flank_range
	deletion_rear_flank = int(row_l[2]) + flank_range
	annotation = []
	gene_dic_out = {}
	exon_dic_out = {}
	UTR_dic_out = {}
	acen_dic_out = {}
	simrep_dic_out = {}
	fr_flank_simprep_dic_out = {}
	fr_flank_repmask_dic_out = {}
	rv_flank_simprep_dic_out = {}
	rv_flank_repmask_dic_out = {}
	repmask_dic_out = {}
	if "chr"+row_l[0] in gene_dict:
		for line in gene_dict["chr"+row_l[0]]:
			gene_name = ""
			status = ""
			if int(line[5]) < int(row_l[1]):
				continue
			elif int(row_l[1]) > int(line[4]) and int(row_l[2]) < int(line[5]):
				status = "within"
				gene_name = line[12] + ";" + line[13] + ";" + line[14] + "," + "within"
				gene_dic_out[gene_name] = 1
			elif int(row_l[1]) <= int(line[4]) and int(line[5]) < int(row_l[2]) < int(line[4]):
				status = "partial"
				gene_name = line[12] + ";" + line[13] + ";" + line[14] + "," + "partial"
				gene_dic_out[gene_name] = 1
			elif int(line[5]) > int(row_l[1]) > int(line[4]) and int(row_l[2]) >= int(line[5]):
				status = "partial"
				gene_name = line[12] + ";" + line[13] + ";" + line[14] + "," + "partial"
				gene_dic_out[gene_name] = 1
			elif int(row_l[1]) <= int(line[4]) and int(row_l[2]) >= int(line[5]):
				status = "cover"
				gene_name = line[12] + ";" + line[13] + ";" + line[14] + "," + "cover"
				gene_dic_out[gene_name] = 1
			if status == "within" or status == "partial" or status == "cover":
				UTR5_list = [line[4], line[6]]
				UTR3_list = [line[7], line[5]]
				UTR = ""
				if int(row_l[1]) > int(UTR5_list[0]) and int(row_l[2]) < int(UTR5_list[1]):
					try:
						utr_prop = str(round((int(row_l[2]) - int(row_l[1]))/(int(UTR5_list[1]) - int(UTR5_list[0])), 4))
					except ZeroDivisionError:
						utr_prop = "non_coding_gene"
					UTR = "5_UTR" + "," + line[12] + "," + utr_prop
					UTR_dic_out[UTR] = 1
				elif int(UTR5_list[0]) >= int(row_l[1]) and int(UTR5_list[1]) > int(row_l[2]) > int(UTR5_list[0]):
					try:
						utr_prop = str(round((int(row_l[2]) - int(row_l[1]))/(int(UTR5_list[1]) - int(UTR5_list[0])), 4))
					except ZeroDivisionError:
						utr_prop = "non_coding_gene"
					utr_prop = str(round((int(row_l[2]) - int(UTR5_list[0]))/(int(UTR5_list[1]) - int(UTR5_list[0])), 4))
					UTR = "5_UTR" + "," + line[12] + "," + utr_prop
					UTR_dic_out[UTR] = 1
				elif int(UTR5_list[0]) < int(row_l[1]) < int(UTR5_list[1]) and int(row_l[2]) >= int(UTR5_list[1]):
					try:
						utr_prop = str(round((int(UTR5_list[1]) - int(row_l[1]))/(int(UTR5_list[1]) - int(UTR5_list[0])), 4))
					except ZeroDivisionError:
						utr_prop = "non_coding_gene"
					UTR = "5_UTR" + "," + line[12] + "," + utr_prop
					UTR_dic_out[UTR] = 1
				elif int(row_l[1]) <= int(UTR5_list[0]) and int(row_l[2]) >= int(UTR5_list[1]):
					try:
						utr_prop = str(round((int(row_l[2]) - int(row_l[1]))/(int(UTR5_list[1]) - int(UTR5_list[0])), 4))
					except ZeroDivisionError:
						utr_prop = "non_coding_gene"
					UTR = "5_UTR" + "," + line[12] + "," + utr_prop
					UTR_dic_out[UTR] = 1
				elif int(row_l[1]) > int(UTR3_list[0]) and int(row_l[2]) < int(UTR3_list[1]):
					try:
						utr_prop = str(round((int(row_l[2]) - int(row_l[1]))/(int(UTR3_list[1]) - int(UTR3_list[0])), 4))
					except ZeroDivisionError:
						utr_prop = "non_coding_gene"
					UTR = "3_UTR" + "," + line[12] + "," + utr_prop
					UTR_dic_out[UTR] = 1
				elif int(UTR3_list[0]) >= int(row_l[1]) and int(UTR3_list[1]) > int(row_l[2]) > int(UTR3_list[0]):
					try:
						utr_prop = str(round((int(row_l[2]) - int(line[7]))/(int(UTR3_list[1]) - int(UTR3_list[0])), 4))
					except ZeroDivisionError:
						utr_prop = "non_coding_gene"
					UTR = "3_UTR" + "," + line[12] + "," + utr_prop
					UTR_dic_out[UTR] = 1
				elif int(UTR3_list[1]) > int(row_l[1]) > int(UTR3_list[0]) and int(row_l[2]) >= int(UTR3_list[1]):
					try:
						utr_prop = str(round((int(UTR3_list[1]) - int(row_l[1]))/(int(UTR3_list[1]) - int(UTR3_list[0])), 4))
					except ZeroDivisionError:
						utr_prop = "non_coding_gene"
					UTR = "3_UTR" + "," + line[12] + "," + utr_prop
					UTR_dic_out[UTR] = 1
				elif int(row_l[1]) <= int(UTR3_list[0]) and int(row_l[2]) >= int(UTR3_list[1]):
					try:
						utr_prop = str(round((int(row_l[2]) - int(row_l[1]))/(int(UTR3_list[1]) - int(UTR3_list[0])), 4))
					except ZeroDivisionError:
						utr_prop = "non_coding_gene"
					UTR = "3_UTR" + "," + line[12] + "," + utr_prop
					UTR_dic_out[UTR] = 1
				exon_stt_list = re.split(",", line[9])
				del exon_stt_list[-1]
				exon_edd_list = re.split(",", line[10])
				del exon_edd_list[-1]
				for stt,edd in zip(exon_stt_list, exon_edd_list):
					exon_pos = ""
					if int(row_l[1]) > int(stt) and int(row_l[2]) < int(edd):
						cover_prop = str(round((int(row_l[2]) - int(row_l[1]))/(int(edd) - int(stt)), 4))
						exon_pos = line[12] + "," + "exon" + "," + stt + "," + edd + "," + "within" + "," + cover_prop
						exon_dic_out[exon_pos] = 1
					elif int(row_l[1]) <= int(stt) and int(edd) > int(row_l[2]) > int(stt):
						cover_prop = str(round((int(row_l[2]) - int(stt))/(int(edd) - int(stt)), 4))
						exon_pos = line[12] + "," + "exon" + "," + stt + "," + edd + "," + "partial" + "," + cover_prop		
						exon_dic_out[exon_pos] = 1
					elif int(edd) > int(row_l[1]) > int(stt) and int(row_l[2]) >= int(edd):
						cover_prop = str(round((int(edd) - int(row_l[1]))/(int(edd) - int(stt)), 4))
						exon_pos = line[12] + "," + "exon" + "," + stt + "," + edd + "," + "partial" + "," + cover_prop
						exon_dic_out[exon_pos] = 1
					elif int(row_l[1]) <= int(stt) and int(row_l[2]) >= int(edd):
						cover_prop = str(round((int(row_l[2]) - int(row_l[1]))/(int(edd) - int(stt)), 4))
						exon_pos = line[12] + "," + "exon" + "," + stt + "," + edd + "," + "cover" + "," + cover_prop
						exon_dic_out[exon_pos] = 1
		if len(gene_dic_out) >= 1:
			gene_annot = ",".join(gene_dic_out.keys())
		else:
			gene_annot = "-"
		if len(exon_dic_out) >= 1:
			exon_annot = ";".join(exon_dic_out.keys())
		else:
			exon_annot = "-"
		if len(UTR_dic_out) >= 1:
			UTR_annot = ",".join(UTR_dic_out.keys())
		else:
			UTR_annot = "-"
	if "chr"+row_l[0] in acen_dict:
		for line in acen_dict["chr"+row_l[0]]:
			acen = ""
			if int(line[1]) < int(row_l[1]):
				continue
			elif int(line[0]) < int(row_l[1]) and int(row_l[2]) < int(line[1]):
				acen_prop = str(round((int(row_l[2]) - int(row_l[1]))/(int(line[1]) - int(line[0])), 4))
				acen = line[2] + "," + "within" + "," + acen_prop
				acen_dic_out[acen] = 1
			elif int(line[0]) >= int(row_l[1]) and int(line[0]) < int(row_l[2]) < int(line[1]):		
				acen_prop = str(round((int(row_l[2]) - int(line[0]))/(int(line[1]) - int(line[0])), 4))
				acen = line[2] + "," + "partial" + "," + acen_prop
				acen_dic_out[acen] = 1
			elif int(line[0]) < int(row_l[1]) < int(line[1]) and int(row_l[2]) >= int(line[1]):
				acen_prop = str(round((int(line[1]) - int(row_l[1]))/(int(line[1]) - int(line[0])), 4))
				acen = line[2] + "," + "partial" + "," + acen_prop
				acen_dic_out[acen] = 1
			elif int(row_l[1]) <= int(line[0]) and int(row_l[2]) >= int(line[1]):
				acen_prop = str(round((int(row_l[2]) - int(row_l[1]))/(int(line[1]) - int(line[0])), 4))
				acen = line[2] + "," + "cover" + "," + acen_prop
				acen_dic_out[acen] = 1
	if len(acen_dic_out) > 0:
		acen_annot = ";".join(acen_dic_out.keys())
	else:
		acen_annot = "-"
	if "chr"+row_l[0] in simplerep_dict:
		for line in simplerep_dict["chr"+row_l[0]]:
			simprep = ""
			fr_flank = ""
			rv_flank = ""
			if int(line[1]) < int(row_l[1]):
				continue
			elif int(line[0]) < int(row_l[1]) and int(row_l[2]) < int(line[1]): 
				simrep_prop = str(round((int(row_l[2]) - int(row_l[1]))/(int(line[1]) - int(line[0])), 4))
				simrep = "simple_repeat" + "," + "within" + "," + simrep_prop + "," + line[2]
				simrep_dic_out[simrep] = 1
			elif int(line[0]) >= int(row_l[1]) and int(line[0]) < int(row_l[2]) < int(line[1]):
				simrep_prop = str(round((int(row_l[2]) - int(line[0]))/(int(row_l[2]) - int(row_l[1])), 4))
				simrep = "simple_repeat" + "," + "partial" + "," + simrep_prop + "," + line[2]
				simrep_dic_out[simrep] = 1
			elif int(line[0]) < int(row_l[1]) < int(line[1]) and int(row_l[2]) >= int(line[1]):
				simrep_prop = str(round((int(line[1]) - int(row_l[1]))/(int(row_l[2]) - int(row_l[1])), 4))
				simrep = "simple_repeat" + "," + "partial" + "," + simrep_prop + "," + line[2]
				simrep_dic_out[simrep] = 1
			elif int(row_l[1]) <= int(line[0]) and int(row_l[2]) >= int(line[1]):
				simrep_prop = str(round((int(row_l[2]) - int(row_l[1]))/(int(line[1]) - int(line[0])), 4))
				simrep = "simple_repeat" + "," + "cover" + "," + simrep_prop + "," + line[2]
				simrep_dic_out[simrep] = 1
			if int(line[1]) < int(row_l[1]) and int(line[0]) > deletion_front_flank:
				fr_flank_prop = str(round((int(line[1]) - int(line[0]))/(int(row_l[1]) - deletion_front_flank), 4))
				fr_flank = "fr_flank" + "," + "covers_simrep" + "," + fr_flank_prop
				fr_flank_simprep_dic_out[fr_flank] = 1
			elif int(line[0]) <= deletion_front_flank and int(row_l[2]) > int(line[1]) >= int(row_l[1]):
				fr_flank_prop = str(round((int(line[1]) - int(line[0]))/(int(row_l[1]) - deletion_front_flank), 4))
				fr_flank = "fr_flank" + "," + "within_simrep" + "," + fr_flank_prop
				fr_flank_simprep_dic_out[fr_flank] = 1
			elif deletion_front_flank < int(line[0]) < int(row_l[1]) and int(row_l[2]) > int(line[1]) >= int(row_l[1]):
				fr_flank_prop = str(round((int(row_l[1]) - int(line[0]))/(int(row_l[1]) - deletion_front_flank), 4))
				fr_flank = "fr_flank" + "," + "partial_simrep" + "," + fr_flank_prop
				fr_flank_simprep_dic_out[fr_flank] = 1
			elif int(line[0]) <= deletion_front_flank and deletion_front_flank < int(line[1]) < int(row_l[1]):
				fr_flank_prop = str(round((int(line[1]) - deletion_front_flank)/(int(row_l[1]) - deletion_front_flank), 4))
				fr_flank = "fr_flank" + "," + "partial_simrep" +  "," + fr_flank_prop
				fr_flank_simprep_dic_out[fr_flank] = 1
			elif int(line[0]) > int(row_l[2]) and int(line[1]) < deletion_rear_flank:
				rv_flank_prop = str(round((int(line[1]) - int(line[0]))/(deletion_rear_flank - int(row_l[2])), 4))
				rv_flank = "rv_flank" + "," + "covers_simrep" + "," + rv_flank_prop
				rv_flank_simprep_dic_out[rv_flank] = 1
			elif int(row_l[1]) < int(line[0]) <= int(row_l[2]) and int(line[1]) >= deletion_rear_flank:
				rv_flank_prop = str(round((int(line[1]) - int(line[0]))/(deletion_rear_flank - int(row_l[2])), 4))
				rv_flank = "rv_flank" + "," + "within_simrep" + "," + rv_flank_prop
				rv_flank_simprep_dic_out[rv_flank] = 1
			elif deletion_rear_flank > int(line[0]) > int(row_l[2]) and int(line[1]) >= deletion_rear_flank:
				rv_flank_prop = str(round((deletion_rear_flank - int(line[0]))/(deletion_rear_flank - int(row_l[2])), 4))
				rv_flank = "rv_flank" + "," + "partial_simrep" + "," + rv_flank_prop
				rv_flank_simprep_dic_out[rv_flank] = 1
			elif int(row_l[1]) < int(line[0]) <= int(row_l[2]) and int(row_l[2]) < int(line[1]) > deletion_rear_flank: 
				rv_flank_prop = str(round((int(line[1]) - int(row_l[2]))/(deletion_rear_flank - int(row_l[2])), 4))
				rv_flank = "rv_flank" + "," + "partial_simrep" + "," + rv_flank_prop
				rv_flank_simprep_dic_out[rv_flank] = 1
	if len(simrep_dic_out) > 0:
		simrep_annot = ";".join(simrep_dic_out.keys())
	else:
		simrep_annot = "-"
	if len(fr_flank_simprep_dic_out) > 0:
		fr_flank_simprep_annot = ";".join(fr_flank_simprep_dic_out.keys())
	else:
		fr_flank_simprep_annot = "-"
	if len(rv_flank_simprep_dic_out) > 0:
		rv_flank_simprep_annot = ";".join(rv_flank_simprep_dic_out.keys())
	else:
		rv_flank_simprep_annot = "-"
	if "chr"+row_l[0] in repmask_dict:
		for line in repmask_dict["chr"+row_l[0]]:
			repmask = "" 
			fr_flank = ""
			rv_flank = ""
			if int(line[1]) < int(row_l[1]):
				continue
			elif int(line[0]) < int(row_l[1]) and int(row_l[2]) < int(line[1]): 
				repmask_prop = str(round((int(row_l[2]) - int(row_l[1]))/(int(line[1]) - int(line[0])), 4))
				repmask = line[4] + "," + line[5] + "," + "within" + "," + repmask_prop
				repmask_dic_out[repmask] = 1
			elif int(line[0]) >= int(row_l[1]) and int(line[0]) < int(row_l[2]) < int(line[1]):
				repmask_prop = str(round((int(row_l[2]) - int(line[0]))/(int(row_l[2]) - int(row_l[1])), 4))
				repmask = line[4] + "," + line[5] + "," + "partial" + "," + repmask_prop
				repmask_dic_out[repmask] = 1
			elif int(line[0]) < int(row_l[1]) < int(line[1]) and int(row_l[2]) >= int(line[1]):
				repmask_prop = str(round((int(line[1]) - int(row_l[1]))/(int(row_l[2]) - int(row_l[1])), 4))
				repmask = line[4] + "," + line[5] + "," + "partial" + "," + repmask_prop
				repmask_dic_out[repmask] = 1
			elif int(row_l[1]) <= int(line[0]) and int(row_l[2]) >= int(line[1]):
				repmask_prop = str(round((int(row_l[2]) - int(row_l[1]))/(int(line[1]) - int(line[0])), 4))
				repmask = line[4] + "," + line[5] + "," + "cover" + "," + repmask_prop
				repmask_dic_out[repmask] = 1
			if int(line[1]) < int(row_l[1]) and int(line[0]) > deletion_front_flank:
				fr_flank_prop = str(round((int(line[1]) - int(line[0]))/(int(row_l[1]) - deletion_front_flank), 4))
				fr_flank = "fr_flank" + "," + "covers_repmask" + "," + line[5] + "," + fr_flank_prop
				fr_flank_repmask_dic_out[fr_flank] = 1
			elif int(line[0]) <= deletion_front_flank and int(row_l[2]) > int(line[1]) >= int(row_l[1]):
				fr_flank_prop = str(round((int(line[1]) - int(line[0]))/(int(row_l[1]) - deletion_front_flank), 4))
				fr_flank = "fr_flank" + "," + "within_repmask" + "," + line[5] + "," + fr_flank_prop
				fr_flank_repmask_dic_out[fr_flank] = 1
			elif deletion_front_flank < int(line[0]) < int(row_l[1]) and int(row_l[2]) > int(line[1]) >= int(row_l[1]):
				fr_flank_prop = str(round((int(row_l[1]) - int(line[0]))/(int(row_l[1]) - deletion_front_flank), 4))
				fr_flank = "fr_flank" + "," + "partial_repmask" + "," + line[5] + "," + fr_flank_prop
				fr_flank_repmask_dic_out[fr_flank] = 1
			elif int(line[0]) <= deletion_front_flank and deletion_front_flank < int(line[1]) < int(row_l[1]):
				fr_flank_prop = str(round((int(line[1]) - deletion_front_flank)/(int(row_l[1]) - deletion_front_flank), 4))
				fr_flank = "fr_flank" + "," + "partial_repmask" + "," + line[5] + "," + fr_flank_prop
				fr_flank_repmask_dic_out[fr_flank] = 1
			elif int(line[0]) > int(row_l[2]) and int(line[1]) < deletion_rear_flank:
				rv_flank_prop = str(round((int(line[1]) - int(line[0]))/(deletion_rear_flank - int(row_l[2])), 4))
				rv_flank = "rv_flank" + "," + "covers_repmask" + "," + line[5] + "," + rv_flank_prop
				rv_flank_repmask_dic_out[rv_flank] = 1
			elif int(row_l[1]) < int(line[0]) <= int(row_l[2]) and int(line[1]) >= deletion_rear_flank:
				rv_flank_prop = str(round((int(line[1]) - int(line[0]))/(deletion_rear_flank - int(row_l[2])), 4))
				rv_flank = "rv_flank" + "," + "within_repmask" + "," + line[5] + "," + rv_flank_prop
				rv_flank_repmask_dic_out[rv_flank] = 1
			elif deletion_rear_flank > int(line[0]) > int(row_l[2]) and int(line[1]) >= deletion_rear_flank:
				rv_flank_prop = str(round((deletion_rear_flank - int(line[0]))/(deletion_rear_flank - int(row_l[2])), 4))
				rv_flank = "rv_flank" + "," + "partial_repmask" + "," + line[5] + "," + rv_flank_prop
				rv_flank_repmask_dic_out[rv_flank] = 1
			elif int(row_l[1]) < int(line[0]) <= int(row_l[2]) and int(row_l[2]) < int(line[1]) < deletion_rear_flank:
				rv_flank_prop = str(round((int(line[1]) - int(row_l[2]))/(deletion_rear_flank - int(row_l[2])), 4))
				rv_flank = "rv_flank" + "," + "partial_repmask" + "," + line[5] + "," + rv_flank_prop
				rv_flank_repmask_dic_out[rv_flank] = 1
	if len(repmask_dic_out) >= 1:
		repmask_annot = ";".join(repmask_dic_out.keys())
	else:
		repmask_annot = "-"
	if len(fr_flank_repmask_dic_out) > 0:
		fr_flank_repmask_annot = ";".join(fr_flank_repmask_dic_out.keys())
	else:
		fr_flank_repmask_annot = "-"
	if len(rv_flank_repmask_dic_out) > 0:
		rv_flank_repmask_annot = ";".join(rv_flank_repmask_dic_out.keys())
	else:
		rv_flank_repmask_annot = "-"
	print("\t".join(row_l[0:3]), gene_annot, exon_annot, UTR_annot, acen_annot, simrep_annot, repmask_annot, fr_flank_simprep_annot, rv_flank_simprep_annot, fr_flank_repmask_annot, rv_flank_repmask_annot, "\t".join(row_l[3:]), sep="\t")
