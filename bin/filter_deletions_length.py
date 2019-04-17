import sys
import re

filtered_file = open(sys.argv[1])
del_len_threshold = int(sys.argv[2])
exclude_len_num = 0
include_num = 0
none_stat = 0
for line in filtered_file:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "chr":
		print("Status", line, sep="\t")
		continue
	del_len = int(line_l[2]) - int(line_l[1])
	len_stat = "NONE"
	samp_num_stat = "NONE"
	if del_len < del_len_threshold:
		len_stat = "Exclude_len"
		exclude_len_num += 1
	elif del_len >= del_len_threshold:
		len_stat = "Include"
		include_num += 1
	else:
		len_stat = "NONE"
		none_stat += 1
	print(len_stat, line, sep="\t")

