import sys
import re

intermediate_file = open(sys.argv[1])
for line in intermediate_file:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "Deletion_status":
		print("\t".join(line_l[1:]))
		continue
	if line_l[0] == "Include":
		print("\t".join(line_l[1:]))

