import sys
import re
import datetime

del_cand_f = open(sys.argv[1])
date = datetime.datetime.today().strftime('%Y%m%d')
reference = sys.argv[2]
assembly = sys.argv[3]
for line in del_cand_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "chr":
		print("##fileformat=VCFv4.1")
		print("##fileDate=" + date)
		print("##reference=" + reference)
		print("##assembly=" + assembly)
		print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">")
		print("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">")
		print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">")
		print("##ALT=<ID=DEL,Description=\"Deletion\">")
		print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
		print("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "\t".join(line_l[29:]), sep="\t") 
		continue
	chrom = line_l[0]
	pos = line_l[1]
	ref = line_l[4]
	end = line_l[2]
	length = line_l[3]
	info_str = "SVTYPE=DEL" + ";" + "END=" + end + ";" + "SVLEN=-" + length
	samples = line_l[31:]
	sample_GT = []
	for sample in samples:
		zygos = ""
		GT = ""
		GT_str = ""
		if sample == "-":
			zygos = "HomoWT"
			GT = "0/0"
			GT_str = GT
		else:
			samp_l = re.split(",", sample)
			if len(samp_l) == 1:
				for sam in samp_l:
					sam_l = re.split(";", sam)
					zygos = sam_l[2]
					if zygos == "Hete":
						GT = "1/0"
					elif zygos == "Homo":
						GT = "1/1"
			elif len(samp_l) > 1:
				sam_l = re.split(";", samp_l[0])
				zygos = sam_l[2]
				if zygos == "Hete":
					GT = "1/0"
				elif zygos == "Homo":
					GT = "1/1"
			GT_str = GT
		sample_GT.append(GT_str)
	print(chrom, pos, ".", ref, "<DEL>", ".", ".", info_str, "GT", "\t".join(sample_GT), sep="\t")
