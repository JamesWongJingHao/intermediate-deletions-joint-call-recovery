import sys
import re
dele_annot_f = open(sys.argv[1])
simrep_prop_coverage_threshold = float(sys.argv[2])
flanking_repeat_coverage_prop_threshold = float(sys.argv[3])

def simple_repeat_eligibility(replist):
	status = "NONE"
	simrep_eligi_list = []
	for rep in replist:
		if "partial" in rep:
			rep_l = re.split(",", rep)
			repeat_prop_simrep = float(rep_l[-2])
			repeat_class = rep_l[0]
			if repeat_class == "simple_repeat" and repeat_prop_simrep >= simrep_prop_coverage_threshold:
				status = "Exclude"
				simrep_eligi_list.append(status)
			else:
				status = "Include"
				simrep_eligi_list.append(status)
	if "Exclude" in simrep_eligi_list:
		status = "Exclude"
	else:
		status = "Include"
	return status
def repeat_masker_eligibility(replist):
	status = "NONE"
	repmask_eligi_list = []
	for rep in replist:
		if "within" in rep:
			rep_l = re.split(",", rep)
			repeat_class = rep_l[0]
			if repeat_class == "Simple_repeat" or repeat_class == "Low_complexity" or repeat_class == "Satellite":
				status = "Exclude"
				repmask_eligi_list.append(status)
			else:
				status = "Include"
				repmask_eligi_list.append(status)
		if "partial" in rep:
			rep_l = re.split(",", rep)
			repeat_prop_repmask = float(rep_l[-1])
			repeat_class = rep_l[0]
			if (repeat_class == "Simple_repeat" or repeat_class == "Low_complexity" or repeat_class == "Satellite") and repeat_prop_repmask >= simrep_prop_coverage_threshold:
				status = "Exclude"
				repmask_eligi_list.append(status)
			else:
				status = "Include"
				repmask_eligi_list.append(status)

	if "Exclude" in repmask_eligi_list:
		status = "Exclude"
	else:
		status = "Include"
	return status
def get_flank_type_list(flank):
	flank_type_list = []
	for type in flank:
		type_l = re.split(",", type)
		flank_type = type_l[1]
		flank_type_list.append(flank_type)
	return flank_type_list
def get_flank_prop_list(flank):
	flank_prop_list = []
	for prop in flank:
		prop_l = re.split(",", prop)
		flank_prop = prop_l[-1]
		flank_prop_list.append(flank_prop)
	return flank_prop_list
def get_flank_partial_props(flank):
	partial_prop_list = []
	for stat in flank:
		if "partial" in stat:
			type_l = re.split(",", stat)
			prop = type_l[-1]
			partial_prop_list.append(prop)
		else:
			a = 1
	return partial_prop_list
def determine_flank_simrep_eligibility(flank_prop_list):
	status = "NONE"
	flank_simrep_eligi_list = []
	for prop in flank_prop_list:
		if float(prop) < flanking_repeat_coverage_prop_threshold:
			status = "Include"
			flank_simrep_eligi_list.append(status)
		if float(prop) >= flanking_repeat_coverage_prop_threshold:
			status = "Exclude"
			flank_simrep_eligi_list.append(status)
	if "Exclude" in flank_simrep_eligi_list:
		status = "Exclude"
	else:
		status = "Include"
	return status 
def determine_flank_repmask_eligibility(repmask_list):
	status = "NONE"
	flank_repmask_eligi_list = []
	for repmask in repmask_list:
		repmask_l = re.split(",", repmask)
		flank_repmask_class = repmask_l[2]
		flank_repmask_type = repmask_l[1]
		flank_repmask_prop = float(repmask_l[-1])
		if flank_repmask_type == "within_repmask" and (flank_repmask_class == "Simple_repeat" or flank_repmask_class == "Low_complexity" or flank_repmask_class == "Satellite"):
			status = "Exclude"
			flank_repmask_eligi_list.append(status)
		elif flank_repmask_type == "partial_repmask" and (flank_repmask_class == "Simple_repeat" or flank_repmask_class == "Low_complexity" or flank_repmask_class == "Satellite"):
			if flank_repmask_prop >= flanking_repeat_coverage_prop_threshold:
				status = "Exclude"
				flank_repmask_eligi_list.append(status)
			else:
				status = "Include"
				flank_repmask_eligi_list.append(status)
		else:
			status = "Include"
			flank_repmask_eligi_list.append(status)
	if "Exclude" in flank_repmask_eligi_list:
		status = "Exclude"
	else:
		status = "Include"
	return status			

for line in dele_annot_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "chr":
		print("Deletion_status", line, sep="\t")
		continue
	dele_length = int(line_l[2]) - int(line_l[1])
	inclusion_status = "NONE"
	final_inclusion_status = []
	sim_rep_stat = line_l[7]
	simple_rep_type_list = []
	simple_rep_prop_list = []
	if not sim_rep_stat == "-":
		sim_rep_l = re.split(";", sim_rep_stat)
		for sim in sim_rep_l:
			rep = re.split(",", sim)
			rep_type = rep[1]
			simple_rep_type_list.append(rep_type)
			rep_prop = rep[-2]
			simple_rep_prop_list.append(rep_prop)
	rep_mask_stat = line_l[8]
	rep_mask_type_list = []
	rep_mask_prop_list = []
	if not rep_mask_stat == "-":
		rep_mask_l = re.split(";", rep_mask_stat)
		for repeat in rep_mask_l:
			rep_mask = re.split(",", repeat)
			rep_mask_type = rep_mask[2]
			rep_mask_type_list.append(rep_mask_type)
			rep_mask_prop = rep_mask[-1]
			rep_mask_prop_list.append(rep_mask_prop)
	fw_flank_simrep_stat = line_l[9]
	fw_flank_simrep_type_list = []
	fw_flank_simrep_prop_list = []
	partial_prop_fw_flank_simrep = []
	if not fw_flank_simrep_stat == "-":
		flank = re.split(";", fw_flank_simrep_stat)
		fw_flank_simrep_type_list = get_flank_type_list(flank)
		fw_flank_simrep_prop_list = get_flank_prop_list(flank)
		partial_prop_fw_flank_simrep = get_flank_partial_props(flank)
	rv_flank_simrep_stat = line_l[10]
	rv_flank_simrep_type_list = []
	rv_flank_simrep_prop_list = []
	partial_prop_rv_flank_simrep = []
	if not rv_flank_simrep_stat == "-":
		flank = re.split(";", rv_flank_simrep_stat)
		rv_flank_simrep_type_list = get_flank_type_list(flank)
		rv_flank_simrep_prop_list = get_flank_prop_list(flank)
		partial_prop_rv_flank_simrep = get_flank_partial_props(flank)
	fw_flank_repmask_stat = line_l[11]
	fw_flank_repmask_type_list = []
	fw_flank_repmask_prop_list = []
	partial_prop_fw_flank_repmask = []
	if not fw_flank_repmask_stat == "-":
		flank = re.split(";", fw_flank_repmask_stat)
		fw_flank_repmask_type_list = get_flank_type_list(flank)
		fw_flank_repmask_prop_list = get_flank_prop_list(flank)
		partial_prop_fw_flank_repmask = get_flank_partial_props(flank)
	rv_flank_repmask_stat = line_l[12]
	rv_flank_repmask_type_list = []
	rv_flank_repmask_prop_list = []
	partial_prop_rv_flank_repmask = []
	if not rv_flank_repmask_stat == "-":
		flank = re.split(";", rv_flank_repmask_stat)
		rv_flank_repmask_type_list = get_flank_type_list(flank)
		rv_flank_repmask_prop_list = get_flank_prop_list(flank)
		partial_prop_rv_flank_repmask = get_flank_partial_props(flank)
	simrep_status = []
	if sim_rep_stat == "-":
		inclusion_status = "Include"
		simrep_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	elif "within" in simple_rep_type_list:
		inclusion_status = "Exclude"
		simrep_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	elif "partial" in simple_rep_type_list:
		replist = re.split(";", sim_rep_stat)
		inclusion_status = simple_repeat_eligibility(replist)
		simrep_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	elif "cover" in simple_rep_type_list:
		inclusion_status = "Include"
		simrep_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	else:
		inclusion_status = "NONE"
		simrep_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	repmask_status = []
	if rep_mask_stat == "-":
		inclusion_status = "Include"
		repmask_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	elif "within" in rep_mask_type_list or "partial" in rep_mask_type_list:
		replist = re.split(";", rep_mask_stat)
		inclusion_status = repeat_masker_eligibility(replist)
		repmask_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	elif "cover" in rep_mask_type_list and not "within" in rep_mask_type_list and not "partial" in rep_mask_type_list:
		inclusion_status = "Include"
		repmask_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	else:
		inclusion_status = "NONE"
		repmask_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	flank_simrep_status = []
	if fw_flank_simrep_stat == "-" and rv_flank_simrep_stat == "-":
		inclusion_status = "Include"
		flank_simrep_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	elif (fw_flank_simrep_stat == "-" and not rv_flank_simrep_stat == "-") or (not fw_flank_simrep_stat == "-" and rv_flank_simrep_stat == "-"):
		inclusion_status = "Include"
		flank_simrep_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	elif "within_simrep" in fw_flank_simrep_type_list and "within_simrep" in rv_flank_simrep_type_list:
		inclusion_status = "Exclude"
		flank_simrep_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	elif "within_simrep" in fw_flank_simrep_type_list and "partial_simrep" in rv_flank_simrep_type_list:
		flank_prop_list = partial_prop_rv_flank_simrep
		inclusion_status = determine_flank_simrep_eligibility(flank_prop_list)
		flank_simrep_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	elif "partial_simrep" in fw_flank_simrep_type_list and "within_simrep" in rv_flank_simrep_type_list:
		flank_prop_list = partial_prop_fw_flank_simrep
		inclusion_status = determine_flank_simrep_eligibility(flank_prop_list)
		flank_simrep_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	elif "partial_simrep" in fw_flank_simrep_type_list and "partial_simrep" in rv_flank_simrep_type_list:
		flank_prop_list = partial_prop_fw_flank_simrep
		inclusion_status = determine_flank_simrep_eligibility(flank_prop_list)
		flank_simrep_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
		flank_prop_list = partial_prop_rv_flank_simrep
		inclusion_status = determine_flank_simrep_eligibility(flank_prop_list)
		flank_simrep_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	elif "covers_simrep" in fw_flank_simrep_type_list or "covers_simrep" in rv_flank_simrep_type_list or ("covers_simrep" in fw_flank_simrep_type_list and "covers_simrep" in rv_flank_simrep_type_list) and (not "within_simrep" in fw_flank_simrep_type_list and not "within_simrep" in rv_flank_simrep_type_list) and not ("partial_simrep" in fw_flank_simrep_type_list and not "partial_simrep" in rv_flank_simrep_type_list):
		inclusion_status = "Include"
		flank_simrep_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	else:
		inclusion_status = "NONE"
		flank_simrep_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	flank_repmask_status = []
	within_repmasks_list = []
	partial_repmask_list = []
	flank_repmask_list = []
	if fw_flank_repmask_stat == "-" and rv_flank_repmask_stat == "-":
		inclusion_status = "Include" 
		flank_repmask_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	elif (fw_flank_repmask_stat == "-" and not rv_flank_repmask_stat == "-") or (not fw_flank_repmask_stat == "-" and rv_flank_repmask_stat == "-"):
		inclusion_status = "Include"
		flank_repmask_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	elif "within_repmask" in fw_flank_repmask_type_list and "within_repmask" in rv_flank_repmask_type_list:
		repmask_list = re.split(";", fw_flank_repmask_stat)
		inclusion_status = determine_flank_repmask_eligibility(repmask_list)
		within_repmasks_list.append(inclusion_status)
		repmask_list = re.split(";", rv_flank_repmask_stat)
		inclusion_status = determine_flank_repmask_eligibility(repmask_list)
		within_repmasks_list.append(inclusion_status)
		if "Include" in within_repmasks_list:
			inclusion_status = "Include"
		elif "Exclude" in within_repmasks_list and not "Include" in within_repmasks_list:
			inclusion_status = "Exclude"
		else:
			inclusion_status = "NONE"
		flank_repmask_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	elif "partial_repmask" in fw_flank_repmask_type_list and "partial_repmask" in rv_flank_repmask_type_list:
		repmask_list = re.split(";", fw_flank_repmask_stat)
		inclusion_status = determine_flank_repmask_eligibility(repmask_list)
		partial_repmask_list.append(inclusion_status)
		repmask_list = re.split(";", rv_flank_repmask_stat)
		inclusion_status = determine_flank_repmask_eligibility(repmask_list)
		partial_repmask_list.append(inclusion_status)
		if "Include" in partial_repmask_list:
			inclusion_status = "Include"
		elif "Exclude" in partial_repmask_list and not "Include" in partial_repmask_list:
			inclusion_status = "Exclude"
		else:
			inclusion_status = "NONE"
		flank_repmask_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	elif ("within_repmask" in fw_flank_repmask_type_list and "partial_repmask" in rv_flank_repmask_type_list) or ("partial_repmask" in fw_flank_repmask_type_list and "within_repmask" in rv_flank_repmask_type_list):
		repmask_list = re.split(";", fw_flank_repmask_stat)
		inclusion_status = determine_flank_repmask_eligibility(repmask_list)
		flank_repmask_list.append(inclusion_status)
		repmask_list = re.split(";", rv_flank_repmask_stat)
		inclusion_status = determine_flank_repmask_eligibility(repmask_list)
		flank_repmask_list.append(inclusion_status)
		if "Include" in flank_repmask_list:
			inclusion_status = "Include"
		elif "Exclude" in flank_repmask_list and not "Include" in flank_repmask_list:
			inclusion_status = "Exclude"
		else:
			inclusion_status = "NONE"
		flank_repmask_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	elif "covers_repmask" in fw_flank_repmask_type_list or "covers_repmask" in rv_flank_repmask_type_list or ("covers_repmask" in fw_flank_repmask_type_list and "covers_repmask" in rv_flank_repmask_type_list) and (not "within_repmask" in fw_flank_repmask_type_list or not "within_repmask" in rv_flank_repmask_type_list) and (not "partial_repmask" in fw_flank_repmask_type_list or not "partial_repmask" in rv_flank_repmask_type_list):
		inclusion_status = "Include"
		flank_repmask_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	else:
		inclusion_status = "NONE"
		flank_repmask_status.append(inclusion_status)
		final_inclusion_status.append(inclusion_status)
	if "Exclude" in final_inclusion_status:
		inclusion_status = "Exclude"
	elif "Include" in final_inclusion_status and not "Exclude" in final_inclusion_status:
		inclusion_status = "Include"
	elif "NONE" in final_inclusion_status:
		inclusion_status = "NONE"
	else:
		inclusion_status = "NONE"
	print(inclusion_status, line, sep="\t")

