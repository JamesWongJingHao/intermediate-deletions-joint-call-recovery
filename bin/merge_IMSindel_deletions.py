import sys
import re
from operator import itemgetter

def define_BP(del_group_list, bp_range_end):
	del_group_list.sort(key=itemgetter(8))	
	grouped_del = {}	
	all_grp_list = []	
	prev = []
	del_group_list2 = []
	for del_tmp in del_group_list:
		if len(prev) and int(del_tmp[8]) - int(prev[8]) > bp_range_end:
			all_grp_list.append(del_group_list2)
			del_group_list2 = []
		del_group_list2.append(del_tmp)
		prev = del_tmp
	if len(del_group_list2):
		all_grp_list.append(del_group_list2)
	indv_output_str = ""
	overall_del_str = ""
	del_cand = {}
	for name in sample_list:
		for grp in all_grp_list:
			overall_del = [grp[0][3], grp[0][4], grp[0][5]]
			overall_del_str = "\t".join(overall_del)
			indv_output_str = "-"
			namesdict = {}
			for indv in grp:
				if name in indv:
					indv_output = [indv[0], indv[1], indv[2], indv[3], indv[4], indv[5], indv[6], indv[8], indv[9], indv[10]]
					indv_output_str = ";".join(indv_output)
					if not name in namesdict:
						namesdict.setdefault(name, []).append(indv_output_str)
					elif name in namesdict:
						namesdict[name].append(indv_output_str)	
				for names in namesdict:
					if len(namesdict[names]) > 1:
						indv_output_str = ",".join(namesdict[names])
					elif len(namesdict[names]) == 1:
						indv_output_str = "".join(namesdict[names])
			if overall_del_str in del_cand.keys():
				del_cand.setdefault(overall_del_str, []).append(indv_output_str)
			elif not overall_del_str in del_cand.keys():
				del_cand[overall_del_str] = [indv_output_str]
	for dele in sorted(del_cand):
		print(dele, "\t".join(del_cand[dele]), sep="\t")
	return del_cand
overall_f = open(sys.argv[1])
bp_range_start = int(sys.argv[2])
bp_range_end = int(sys.argv[3])
sample_list = []
for line in overall_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "sample":
		continue
	if not line_l[0] in sample_list:
		sample_list.append(line_l[0])
sample_list = sorted(sample_list) 		
print("chr", "stt", "edd", "\t".join(sample_list), sep="\t")
overall_f = open(sys.argv[1])
bp_range_start = int(sys.argv[2])
bp_range_end = int(sys.argv[3])
previous = []
del_group_list = []
all_grp_list = []
del_cand = {}
for line in overall_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "sample":
		continue
	if len(previous) > 0:
		if line_l[4] != previous[4] or (line_l[4] == previous[4] and abs(int(line_l[5]) - int(previous[5])) > bp_range_start):
			if len(del_group_list) >= 2:
				all_grp_list = define_BP(del_group_list, bp_range_end)
			else:
				del_group_list = define_BP(del_group_list, bp_range_end)
			del_group_list = []
	del_group_list.append(line_l)
	previous = line_l	
all_grp_list = define_BP(del_group_list, bp_range_end)

