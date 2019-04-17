import sys
import re
import pysam

HWE_filt_f = open(sys.argv[1])
bam_dir_f = open(sys.argv[2])
extract_range = int(sys.argv[3])
bp_basequal_extract_region = int(sys.argv[4])
bam_dir_dict = {}
for row in bam_dir_f:
	row = row.replace('\n', '')
	row_l = re.split("\t", row)
	samplename = row_l[0] 
#	contig = row_l[1]
	bamfile_location = row_l[-1]
#	samplename_contig = samplename + "_" + contig
#	if not samplename_contig in bam_dir_dict:
#		bam_dir_dict.setdefault(samplename_contig, []).append)(bamfile_location)
	if not samplename in bam_dir_dict:
		bam_dir_dict.setdefault(samplename, []).append(bamfile_location)
for line in HWE_filt_f:
	line = line.replace('\n', '')
	line_l = re.split("\t", line)
	if line_l[0] == "chr":
		print("\t".join(line_l[0:21]), "Avg_forward_region_readdepth", "Avg_deletion_region_readdepth", "Avg_reverse_region_readdepth", "Avg_forward_reads_softclip_score", "Avg_reverse_reads_softclip_score", "\t".join(line_l[21:]), sep="\t")
		continue
	samples = line_l[23:]
	samples_list = []
	for sample in samples:
		if not sample == "-":
			samples_list.append(sample)
	deletion_name = line_l[0] + "_" + line_l[1] + "_" + line_l[2]
	chrom = line_l[0]
#	deletion_contig = line_l[4]
	del_stt = line_l[1]
	del_edd = line_l[2]
	del_extract_stt = int(del_stt) - extract_range
	del_extract_edd = int(del_edd) + extract_range
	del_region_size = int(del_edd) - int(del_stt)
#	del_extract_stt = del_contig_stt - extract_range
#	del_extract_edd = del_contig_edd + extract_range
#	del_region_size = del_contig_edd - del_contig_stt
	fw_bp_stt = int(del_stt) - del_region_size
	rv_bp_edd = int(del_edd) + del_region_size
	overall_fw_region_depth = 0
	overall_deletion_region_depth = 0
	overall_rv_region_depth = 0
	fw_reads_overall_score = 0
	rv_reads_overall_score = 0
	for sample in samples_list:
		sample_l = re.split(";", sample)
		samplename = sample_l[0]
#		sample_detail = samplename + "_" + deletion_contig
#		sample_bam_file = "".join(bam_dir_dict[sample_detail])
		sample_bam_file = "".join(bam_dir_dict[samplename])
		sample_alignment_file = pysam.AlignmentFile(sample_bam_file)
		fw_bp_region_depth = 0
		avg_fw_bp_region_depth = 0
		for pileupcolumn in sample_alignment_file.pileup(chrom, fw_bp_stt, int(del_stt)):
			fw_bp_region_depth += pileupcolumn.n
			region_size = int(del_stt) - fw_bp_stt
			avg_fw_bp_region_depth = round((fw_bp_region_depth/region_size), 0)
		del_region_depth = 0
		avg_del_region_depth = 0
		for pileupcolumn in sample_alignment_file.pileup(chrom, int(del_stt), int(del_edd)):
			del_region_depth += pileupcolumn.n
			region_size = int(del_edd) - int(del_stt)
			avg_del_region_depth = round((del_region_depth/region_size), 0)
		rv_bp_region_depth = 0
		avg_rv_bp_region_depth = 0
		for pileupcolumn in sample_alignment_file.pileup(chrom, int(del_edd), rv_bp_edd):
			rv_bp_region_depth += pileupcolumn.n
			region_size = rv_bp_edd - int(del_edd)
			avg_rv_bp_region_depth = round((rv_bp_region_depth/region_size), 0)
		sample_region = samplename + "_" + deletion_name
		overall_fw_region_depth += avg_fw_bp_region_depth
		overall_deletion_region_depth += avg_del_region_depth
		overall_rv_region_depth += avg_rv_bp_region_depth
		basequal_stt_bp = int(del_stt) - bp_basequal_extract_region
		basequal_edd_bp = int(del_edd) + bp_basequal_extract_region
		avg_fw_read_S_qscore = 0
		avg_rv_read_S_qscore = 0
		for read in sample_alignment_file.fetch(chrom, basequal_stt_bp, int(del_stt)):
			samread = str(read)
			samread = samread.replace('\n', '')
			samread_l = re.split("\t", samread)
			cigarline = samread_l[5]
			readflag = int(samread_l[1])
			if readflag & 0x400 or readflag & 0x800 or readflag & 0x4:
				continue
			cigar_len = re.split('\D', cigarline)
			cigar_type = re.split('[0-9]+', cigarline)
			if (len(cigar_type) == 2 and cigar_type[-1] == "M") or cigar_type[-1] == "M":
				continue
			if "I" in cigar_type or "D" in cigar_type:
				continue
			unaligned_pos1 = 0	
			if cigar_type[-1] == "S":
				unaligned_pos1 = int(cigar_len[-2])
				qscore_line = samread_l[10]
				score_list = qscore_line.split("[")
				score = score_list[1]
				score_l = re.split("]", score)
				scores = score_l[0:-1]
				scores_l = re.split(",", scores[0])
				scores_a = [int(x) for x in scores_l]
				fw_S_qscore_sum = 0 
				if unaligned_pos1 > 0:
					unaligned1_qscore = scores_a[len(scores_a) - unaligned_pos1:len(scores_a)]
					for qscore in unaligned1_qscore:
						fw_S_qscore_sum += qscore
				avg_fw_read_S_qscore = round((fw_S_qscore_sum/len(unaligned1_qscore)), 0)
		for read in sample_alignment_file.fetch(chrom, int(del_edd), basequal_edd_bp):
			samread = str(read)
			samread = samread.replace('\n', '')
			samread_l = re.split("\t", samread)
			cigarline = samread_l[5]
			readflag = int(samread_l[1])
			if readflag & 0x400 or readflag & 0x800 or readflag & 0x4:
				continue
			cigar_len = re.split('\D', cigarline)
			cigar_type = re.split('[0-9]+', cigarline)
			if (len(cigar_type) == 2 and cigar_type[-1] == "M") or cigar_type[1] == "M":
				continue
			if "I" in cigar_type or "D" in cigar_type:
				continue
			unaligned_pos2 = 0
			if cigar_type[1] == "S":
				unaligned_pos2 = int(cigar_len[0])
				qscore_line = samread_l[10]
				score_list = qscore_line.split("[")
				score = score_list[1]
				score_l = re.split("]", score)
				scores = score_l[0:-1]
				scores_l = re.split(",", scores[0])
				scores_a = [int(x) for x in scores_l]
				rv_S_qscores_sum = 0
				if unaligned_pos2 > 0:
					unaligned2_qscore = scores_a[0:unaligned_pos2]
					for qscore in unaligned2_qscore:
						rv_S_qscores_sum += qscore
				avg_rv_read_S_qscore = round((rv_S_qscores_sum/len(unaligned2_qscore)), 0)
		fw_reads_overall_score += avg_fw_read_S_qscore		
		rv_reads_overall_score += avg_rv_read_S_qscore
	avg_fw_region_depth = round((overall_fw_region_depth/len(samples_list)), 0)
	avg_deletion_region_depth = round((overall_deletion_region_depth/len(samples_list)), 0)
	avg_rv_region_depth = round((overall_rv_region_depth/len(samples_list)), 0)
	avg_fw_region_readscores = round((fw_reads_overall_score/len(samples_list)), 0)
	avg_rv_region_readscores = round((rv_reads_overall_score/len(samples_list)), 0)		
	print("\t".join(line_l[0:23]), avg_fw_region_depth, avg_deletion_region_depth, avg_rv_region_depth, avg_fw_region_readscores, avg_rv_region_readscores, "\t".join(line_l[23:]), sep="\t")



