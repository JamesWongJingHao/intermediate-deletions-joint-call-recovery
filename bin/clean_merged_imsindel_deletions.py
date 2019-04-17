import sys
import re
res_f = open(sys.argv[1])
print("sample", ">indel_type", "call_type", "chr", "sttpos", "endpos", "indel_length", "indel_str", "#indel_depth", "#ttl_depth", "details(indelcall_indeltype_depth)", "clip_sttpos", "depth(>=10)", sep = "\t")
for line in res_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "sample":
		continue
	print(line)


