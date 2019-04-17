import sys
import re

intermdiate_f = open(sys.argv[1])
for line in intermdiate_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "status":
		print("\t".join(line_l[1:]))
		continue
	if not line_l[0] == "Tier_1":
		continue
	print("\t".join(line_l[1:]))		
