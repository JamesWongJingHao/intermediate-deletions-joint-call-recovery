#>indel_type      call_type       chr    contig   sttpos          endpos          indel_length    indel_str       #indel_depth    #ttl_depth      details(indelcall_indeltype_depth)      clip_sttpos     depth(>=10)
#DEL     Hete    gi|13626247|ref|NT_025975.2|    NT_025975.2     10566   14153   3588    TTTTTGCATTCCATTACATTCTATGACATTCGATTCCGTTTCATTGCATTCCATTCCATACATTTTTATTCCATTCGAGACCGTAGCATTCCACTTTATTCCAGGCCTGTCCATTACACTACATTCCCTTCCATTCCAATGAATTCCATTCCATTCCAATCCATTCCTTTCCTTTCGCTTGCATTCCATTCTATTCTCTTCTACTGCATACAATTTCACTCCATTCGTTCCCATTCCATTCAATTCCATTCCA
#DEL     Hete    gi|13626247|ref|NT_025975.2|    NT_025975.2     10732   10732   1       C       7       58      DEL_SID_7_14320 14320   "High,Lowfreq"

#contig.md:
#9606    1       1       1       +       start   na      CONTIG  Primary_Assembly        -2
#9606    1       10001   267719  +       NT_077402.2     224514618       CONTIG  Primary_Assembly        1
#9606    1       317720  471368  +       NT_077912.1     29794409        CONTIG  Primary_Assembly        1
#>indel_type      call_type       chr     sttpos          endpos          indel_length    indel_str       #indel_depth    #ttl_depth      details(indelcall_indeltype_depth)      clip_sttpos     depth(>=10)

#>indel_type      call_type       contig          sttpos          endpos          chr     chrm_stt        chrm_end        indel_length    indel_str       #indel_depth    #ttl_depth      details(indelcall_indeltype_depth)      clip_sttpos     depth(>=10)
#DEL     Hete    gi|13626247|ref|NT_025975.2|    7997    8011    Y       58827358        58827372        15      TCCATTCAATTCCTT 6       435     DEL_SID_6_8002  8002    High,Lowfreq
#DEL     Hete    gi|13626247|ref|NT_025975.2|    8022    8035    Y       58827383        58827396        14      TTGATTTGATTCCA  2       257     DEL_SID_2_8039  8039    High,Lowfreq

import sys
import re
results_f = open(sys.argv[1])
contig_f = open(sys.argv[2])
print(">indel_type", "\t", "call_type", "\t", "contig", "\t", "sttpos", "\t", "endpos", "\t", "chr", "\t", "chrm_stt", "\t", "chrm_end", "\t", "indel_length", "\t", "indel_str", "\t", "#indel_depth", "\t", "#ttl_depth", "\t", "details(indelcall_indeltype_depth)", "\t", "clip_sttpos", "\t", "depth(>=10)")
contig_pos = []
contig_dict = {}
for row in contig_f:
	row = row.replace('\n', '')
	row_l = re.split("\t", row)
	contig_pos = [row_l[1], row_l[2], row_l[3]]
	contig_pos_str = "\t".join(contig_pos)
	contig_dict.setdefault(row_l[5], []).append(contig_pos_str)
contig_strt = 0
chrmstart = 0
chrmend = 0
chrom = ""
for line in results_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == ">indel_type ":
		continue
	results_contig = re.split("\|", line_l[2])
	res_con = results_contig[3]
	for pos in contig_dict[res_con]:
		pos_l = re.split("\t", pos)
		chrom = pos_l[0]
		contig_strt = int(pos_l[1])		
	chrmstart = (contig_strt + int(line_l[3])) - 1
	chrmend = (contig_strt + int(line_l[4])) - 1
	print(line_l[0], line_l[1], line_l[2], line_l[3], line_l[4], chrom, chrmstart, chrmend, line_l[5], line_l[6], line_l[7], line_l[8], line_l[9], line_l[10], line_l[11], sep="\t")
