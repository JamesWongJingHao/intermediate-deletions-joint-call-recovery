import sys
import re

intermediate_file = open(sys.argv[1])
for line in intermediate_file:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "Status":
		print("\t".join(line_l[1:4]), "deletion_length", "\t".join(line_l[4:]), sep="\t")
		continue
	dele_length = int(line_l[3]) - int(line_l[2])
	status = line_l[0]
	if status == "Include":
		print("\t".join(line_l[1:4]), dele_length, "\t".join(line_l[4:]), sep="\t")

