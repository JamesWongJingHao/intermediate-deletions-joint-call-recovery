import sys
import re

deletion_HWE_res_f = open(sys.argv[1])
HWE_exclusion_threshold = float(sys.argv[2])

for line in deletion_HWE_res_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "chr":
		print(line)
		continue
	HWE_res = float(line_l[22])
	if HWE_res < HWE_exclusion_threshold:
		continue
	print(line)

