import sys
import re

res_f = open(sys.argv[1])
print(">indel_type", "\t", "call_type", "\t", "chr", "\t", "sttpos", "\t", "endpos", "\t", "indel_length", "\t", "indel_str", "\t", "#indel_depth", "\t", "#ttl_depth", "\t", "details(indelcall_indeltype_depth)", "\t", "clip_sttpos", "\t", "depth(>=10)")
for line in res_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == ">indel_type":
		continue
	print(line)


