import sys
import re

ref_bases_f = open(sys.argv[1])
del_f = open(sys.argv[2])
pos_ref_dict = {}
for row in ref_bases_f:
	row = row.replace('\n', '')
	row_l = re.split("\t", row)
	if row_l[0] == "sample":
		continue
	del_position = row_l[3] + "_" + row_l[4] + "_" + row_l[5]
	ref_bases = row_l[7]
	if not del_position in pos_ref_dict:
		pos_ref_dict[del_position] = [ref_bases]
	elif del_position in pos_ref_dict and not ref_bases in pos_ref_dict[del_position]:
		pos_ref_dict[del_position].append(ref_bases)
for line in del_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "chr":
		print("\t".join(line_l[0:4]), "ref", "\t".join(line_l[4:]), sep="\t")
		continue
	deletion_pos = line_l[0] + "_" + line_l[1] + "_" + line_l[2]	
	if deletion_pos in pos_ref_dict:
		ref = "".join(pos_ref_dict[deletion_pos])
	else:
		ref = "NA"
	print("\t".join(line_l[0:4]), ref, "\t".join(line_l[4:]), sep="\t")
