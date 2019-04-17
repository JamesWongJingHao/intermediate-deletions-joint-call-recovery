import sys
import re
samplename_f = open(sys.argv[1])
for line in samplename_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "sample":
		print(line)
		continue
	if line_l[1] == "DEL":
		print(line)
