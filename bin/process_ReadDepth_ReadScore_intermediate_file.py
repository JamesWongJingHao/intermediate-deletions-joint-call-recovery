import sys
import re

intermediate_f = open(sys.argv[1])
for line in intermediate_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "Status":
		print("\t".join(line_l[1:]))
		continue
	if not line_l[0] == "Include":
		continue
	print("\t".join(line_l[1:]))
