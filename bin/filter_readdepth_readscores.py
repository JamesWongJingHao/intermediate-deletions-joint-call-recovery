import sys
import re

HWE_RD_RS_f = open(sys.argv[1])
coverage_depth_threshold = int(sys.argv[2])
readscore_threshold = int(sys.argv[3])
for line in HWE_RD_RS_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "chr":
		print("Status", line, sep="\t")
		continue
	fw_RD = int(float(line_l[23]))
	delreg_RD = int(float(line_l[24]))
	rv_RD = int(float(line_l[25]))
	fw_RS = int(float(line_l[26]))
	rv_RS = int(float(line_l[27]))
	status = "UNKNOWN"
	if fw_RD > coverage_depth_threshold and delreg_RD > coverage_depth_threshold and rv_RD > coverage_depth_threshold:
		status = "exclude_overdepth"
	elif delreg_RD > fw_RD and delreg_RD > rv_RD:
		status = "exclude_deletion_overdepth"
	elif fw_RS < readscore_threshold or rv_RS < readscore_threshold:
		status = "exclude_readscore"
	else:
		status = "Include"
	print(status, line, sep="\t")

