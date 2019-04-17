import sys
import re

repmask_f = open(sys.argv[1])
overall_res_f = open(sys.argv[2])
repmask_dict = {}
for line in repmask_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "genoName":
		continue
	chrom = line_l[0]
	rep = [line_l[1], line_l[2], line_l[7], line_l[6], line_l[5]]
	repmask_dict.setdefault(chrom, []).append(rep)
output_str = ""
for row in overall_res_f:
	row = row.replace('\n', '')
	row_l = re.split("\t", row)
	if row_l[0] == "chr":
		print("\t".join(row_l[0:26]), "del_front_rep", "del_rear_rep", "\t".join(row_l[26:]), sep="\t")
		continue
	deletion_stt = int(row_l[1])
	deletion_edd = int(row_l[2])
	annotation = []
	fr_flank_repmask_dic_out = {}
	rv_flank_repmask_dic_out = {}
	repmask_dic_out = {}
	if "chr"+row_l[0] in repmask_dict:
		for line in repmask_dict["chr"+row_l[0]]:
			repmask = "" 
			fr_flank = ""
			rv_flank = ""
			if int(line[1]) < int(row_l[1]):
				continue
			if int(line[0]) <= deletion_stt <= int(line[1]):
				fr_flank = "front" + "," + "within" + "," + line[3] + "," + line[-1]
				fr_flank_repmask_dic_out[fr_flank] = 1
			if int(line[0]) <= deletion_edd <= int(line[1]):
				rv_flank = "rear" + "," + "within" + "," + line[3] + "," + line[-1]
				rv_flank_repmask_dic_out[rv_flank] = 1
	if len(fr_flank_repmask_dic_out) > 0:
		fr_flank_repmask_annot = ";".join(fr_flank_repmask_dic_out.keys())
	else:
		fr_flank_repmask_annot = "-"

	if len(rv_flank_repmask_dic_out) > 0:
		rv_flank_repmask_annot = ";".join(rv_flank_repmask_dic_out.keys())
	else:
		rv_flank_repmask_annot = "-"
	print("\t".join(row_l[0:28]), fr_flank_repmask_annot, rv_flank_repmask_annot, "\t".join(row_l[28:]), sep="\t")
	
