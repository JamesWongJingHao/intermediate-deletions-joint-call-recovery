import sys
import re

input_del_f = open(sys.argv[1])
samplename = sys.argv[2]
for line in input_del_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == ">indel_type ":
		print("sample", line, sep = "\t")
		continue
	sample = samplename
	print(sample, line, sep = "\t")
