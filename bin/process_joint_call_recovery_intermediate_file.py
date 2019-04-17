import sys
import re

inter_f = open(sys.argv[1])
for line in inter_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "Deletion_status":
		print("\t".join(line_l[1:]))
	if not line_l[0] == "Rescue":
		continue
	print("\t".join(line_l[1:]))
