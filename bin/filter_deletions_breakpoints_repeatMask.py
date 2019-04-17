import sys
import re

HWE_filt_f = open(sys.argv[1])
def get_breakpoint_type(bp_stat):
	for bp_stat in bp_stat:
		bp_stat_a = re.split(";", bp_stat)
	for stat in bp_stat_a:
		stat_l = re.split(",", stat)
		rep_type = stat_l[2]
	return rep_type
def get_breakpoint_name(bp_stat):
	for bp_stat in bp_stat:
		bp_stat_a = re.split(";", bp_stat)
	for stat in bp_stat_a:
		stat_l = re.split(",", stat)
		rep_name = stat_l[3]
	return rep_name
def get_deletion_status(rep_type_fw, rep_type_rv, rep_name_fw, rep_name_rv):
	if rep_type_fw == "Simple_repeat" and rep_type_rv == "Simple_repeat" and rep_name_fw == rep_name_rv:
		status = "Tier_1"
	elif (rep_type_fw == rep_type_rv and rep_name_fw == rep_name_rv) and (not "within" in rep_name_stat_dict[rep_type_fw]):
		status = "Tier_2"
	elif (rep_type_fw == rep_type_rv) and ("Alu" in rep_name_fw and "Alu" in rep_name_rv) and (not "within" in rep_name_stat_dict[rep_type_fw]):
		status = "Tier_2"	
	else:
		status = "Tier_1"
	return status
for line in HWE_filt_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "chr":
		print("status", line, sep="\t")
		continue
	fw_bp_stat = line_l[28]
	rv_bp_stat = line_l[29]
	del_repmask_annot = line_l[9]
	status = "NA"
	if (fw_bp_stat == "-" and rv_bp_stat == "-") or (fw_bp_stat == "-" and not rv_bp_stat == "-") or (not fw_bp_stat == "-" and rv_bp_stat == "-"):
		status = "Tier_1"
	elif not fw_bp_stat == "-" and not rv_bp_stat == "-":
		del_repmask = re.split(";", del_repmask_annot)
		rep_name_stat_dict = {}
		for rep in del_repmask:
			rep_l = re.split(",", rep)
			rep_stat = rep_l[2]
			rep_name = rep_l[0]
			rep_class = rep_l[1]
			if not rep_name in rep_name_stat_dict:
				rep_name_stat_dict[rep_name] = [rep_stat]
			else:
				rep_name_stat_dict[rep_name].append(rep_stat)
		fw_bp_stat_a = re.split(";", fw_bp_stat)
		rv_bp_stat_a = re.split(";", rv_bp_stat)
		status_del = []
		if len(fw_bp_stat_a) == 1 and len(rv_bp_stat_a) == 1:
			bp_stat = fw_bp_stat_a
			rep_type_fw = get_breakpoint_type(bp_stat)
			rep_name_fw = get_breakpoint_name(bp_stat)
			bp_stat = rv_bp_stat_a
			rep_type_rv = get_breakpoint_type(bp_stat)
			rep_name_rv = get_breakpoint_name(bp_stat)
			status = get_deletion_status(rep_type_fw, rep_type_rv, rep_name_fw, rep_name_rv)
			status_del.append(status)
		elif len(fw_bp_stat_a) == 1 and len(rv_bp_stat_a) > 1:
			bp_stat = fw_bp_stat_a
			rep_type_fw = get_breakpoint_type(bp_stat)
			rep_name_fw = get_breakpoint_name(bp_stat)
			for stat in rv_bp_stat_a:
				bp_stat = [stat]
				rep_type_rv = get_breakpoint_type(bp_stat)
				rep_name_rv = get_breakpoint_name(bp_stat)
				status = get_deletion_status(rep_type_fw, rep_type_rv, rep_name_fw, rep_name_rv)
				status_del.append(status)
		elif len(fw_bp_stat_a) > 1 and len(rv_bp_stat_a) == 1:
			bp_stat = rv_bp_stat_a
			rep_type_rv = get_breakpoint_type(bp_stat)
			rep_name_rv = get_breakpoint_name(bp_stat)
			for stat in fw_bp_stat_a:
				bp_stat = [stat]
				rep_type_fw = get_breakpoint_type(bp_stat)
				rep_name_fw = get_breakpoint_name(bp_stat)
				status = get_deletion_status(rep_type_fw, rep_type_rv, rep_name_fw, rep_name_rv)
				status_del.append(status)
		if len(status_del) > 1 and "Tier_2" in status_del:
			status = "Tier_2"
	print(status, line, sep="\t")
