#!/bin/bash
#!/usr/bin/r
path=$0
if [ "$path" == "${path%/*}" ]; then src=.; else src=${path%/*}; fi
input_dir=$1
out_dir=$2
date=`date "+%Y%m%d"`
gene_file=$src/required_files/hg19_genes_ucsc.txt
acen_file=$src/required_files/hg19_telomere_centromere.txt
simplerep_file=$src/required_files/hg19_simple_repeats_ucsc.txt
repmask_file=$src/required_files/hg19_repeatmasker_ucsc.txt
bamfiles_list=$src/required_files/bam_files_locations.txt
repmask_BED_file=$src/required_files/hg19_repeatmasker_ucsc_BED.txt
merge_deletion_breakpoints_start=`python $src/bin/read_config_file.py $src/config merge_deletion_breakpoints_start`
merge_deletion_breakpoints_end=`python $src/bin/read_config_file.py $src/config merge_deletion_breakpoints_end`
annotation_flank=`python $src/bin/read_config_file.py $src/config annotation_flank`
simplerepeat_proportion_threshold=`python $src/bin/read_config_file.py $src/config simplerepeat_proportion_threshold`
flanking_repeat_coverage_proportion_threshold=`python $src/bin/read_config_file.py $src/config flanking_repeat_coverage_proportion_threshold`
LD_dels_support_reads_threshold=`python $src/bin/read_config_file.py $src/config LD_dels_support_reads_threshold`
SID_dels_number_threshold=`python $src/bin/read_config_file.py $src/config SID_dels_number_threshold`
LD_dels_number_threshold=`python $src/bin/read_config_file.py $src/config LD_dels_number_threshold`
LD_sr_number_threshold=`python $src/bin/read_config_file.py $src/config LD_sr_number_threshold`
del_len_recovery_threshold=`python $src/bin/read_config_file.py $src/config del_len_recovery_threshold`
del_len_filter_threshold=`python $src/bin/read_config_file.py $src/config del_len_filter_threshold`
HWE_exclusion_threshold=`python $src/bin/read_config_file.py $src/config HWE_exclusion_threshold`
bamfile_extract_range=`python $src/bin/read_config_file.py $src/config bamfile_extract_range`
breakpoint_score_extract_range=`python $src/bin/read_config_file.py $src/config breakpoint_score_extract_range`
readdepth_filter_threshold=`python $src/bin/read_config_file.py $src/config readdepth_filter_threshold`
readscore_filter_threshold=`python $src/bin/read_config_file.py $src/config readscore_filter_threshold`
reference=`python $src/bin/read_config_file.py $src/config reference`
assembly=`python $src/bin/read_config_file.py $src/config assembly`
if [ ! -d $out_dir ]; then mkdir $out_dir; fi

echo "Concatenating deletion lists from individual samples... "
cat $input_dir/*.deletions.* > $out_dir/temporary_concatenated_deletions.$date
python $src/bin/clean_merged_imsindel_deletions.py $out_dir/temporary_concatenated_deletions.$date > $out_dir/temp_joined_deletions_clean.$date
sort -k4,4V -k5,5n $out_dir/temp_joined_deletions_clean.$date > $out_dir/IMSindel_overall_deletions.$date
exit_value=$?; if [ $exit_value != 0 ]; then echo "<< Error >>"; exit 1; else echo "<< Finished >>"; fi
rm $out_dir/temporary_concatenated_deletions.$date
rm $out_dir/temp_joined_deletions_clean.$date

echo "Merging deletion calls... "
python $src/bin/merge_IMSindel_deletions.py $out_dir/IMSindel_overall_deletions.$date $merge_deletion_breakpoints_start $merge_deletion_breakpoints_end > $out_dir/IMSindel.overall.deletions.merged.$date
#rm $out_dir/IMSindel_overall_deletions.$date
exit_value=$?; if [ $exit_value != 0 ]; then echo "<< Error >>"; exit 1; else echo "<< Finished >>"; fi

echo "Annotating deletions... "
python $src/bin/IMSindel_merged_deletions_annotation.py $gene_file $acen_file $simplerep_file $repmask_file $out_dir/IMSindel.overall.deletions.merged.$date $annotation_flank > $out_dir/IMSindel.deletions.annotated.$date
exit_value=$?; if [ $exit_value != 0 ]; then echo "<< Error >>"; exit 1; else echo "<< Finished >>"; fi

echo "Filtering deletions and conducting joint call recovery... "
python $src/bin/annotation_filter_IMSindel_deletions.py $out_dir/IMSindel.deletions.annotated.$date $simplerepeat_proportion_threshold $flanking_repeat_coverage_proportion_threshold > $out_dir/IMSindel.deletions.annotated.filtered.intermediate
python $src/bin/process_annotation_filtered_intermediate_file.py $out_dir/IMSindel.deletions.annotated.filtered.intermediate > $out_dir/IMSindel.deletions.annotated.filtered.$date
python $src/bin/joint_call_recovery.py $out_dir/IMSindel.deletions.annotated.filtered.$date $LD_dels_support_reads_threshold $SID_dels_number_threshold $LD_dels_number_threshold $LD_sr_number_threshold $del_len_recovery_threshold > $out_dir/IMSindel.deletions.annotated.filtered.jointrecov.intermediate
python $src/bin/process_joint_call_recovery_intermediate_file.py $out_dir/IMSindel.deletions.annotated.filtered.jointrecov.intermediate > $out_dir/IMSindel.deletions.annotated.filtered.joint_recovery.$date
python $src/bin/filter_deletions_length.py $out_dir/IMSindel.deletions.annotated.filtered.joint_recovery.$date $del_len_filter_threshold > $out_dir/IMSindel.deletions.annotated.filtered.joint_recovery.length_filtered.intermediate
python $src/bin/process_deletions_length_filter.py $out_dir/IMSindel.deletions.annotated.filtered.joint_recovery.length_filtered.intermediate > $out_dir/IMSindel.deletions.annotated.joint_recovery.filtered.$date
python $src/bin/convert_deletions_dataframe_HWE.py $out_dir/IMSindel.deletions.annotated.joint_recovery.filtered.$date > $out_dir/IMSindel.deletions.dataframe.HWE.txt
R --slave --vanilla $out_dir/IMSindel.deletions.dataframe.HWE.txt < $src/bin/calculate_HWE.R
python $src/bin/merge_deletions_HWE.py $out_dir/IMSindel_deletions_HWE_results.txt $out_dir/IMSindel.deletions.annotated.joint_recovery.filtered.$date > $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE
python $src/bin/filter_by_HWE.py $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE $HWE_exclusion_threshold > $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.$date
python $src/bin/calculate_deletions_ReadDepth_ReadScores_v2.py $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.$date $bamfiles_list $bamfile_extract_range $breakpoint_score_extract_range > $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD_RS.$date
python $src/bin/filter_readdepth_readscores.py $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD_RS.$date $readdepth_filter_threshold $readscore_filter_threshold > $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD.RS.intermediate
python $src/bin/process_ReadDepth_ReadScore_intermediate_file.py $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD.RS.intermediate > $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD.RS.temp 
python $src/bin/annotate_flanks_repeatMask.py $repmask_BED_file $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD.RS.temp > $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD_RS_filtered.$date
exit_value=$?; if [ $exit_value != 0 ]; then echo "<< Error >>"; exit 1; else echo "<< Finished >>"; fi
rm $out_dir/IMSindel.deletions.annotated.filtered.intermediate
rm $out_dir/IMSindel.deletions.annotated.filtered.jointrecov.intermediate 
rm $out_dir/IMSindel.deletions.annotated.filtered.joint_recovery.length_filtered.intermediate
rm $out_dir/IMSindel.deletions.dataframe.HWE.txt
rm $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE
rm $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD.RS.intermediate
rm $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD.RS.temp

echo "Generating high-confidence deletion call set and VCF file... "
python $src/bin/filter_deletions_breakpoints_repeatMask.py $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD_RS_filtered.$date > $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD_RS_filtered.BP.intermediate
python $src/bin/process_deletions_breakpoints_repeatMask_intermediate_file.py $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD_RS_filtered.BP.intermediate > $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD_RS_filtered.BP.intermediate.temp 
python $src/bin/append_deletions_reference_bases.py $out_dir/IMSindel_overall_deletions.$date $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD_RS_filtered.BP.intermediate.temp > $out_dir/intermediate_deletions_HC.$date
python $src/bin/convert_deletions_to_vcf.py $out_dir/intermediate_deletions_HC.$date $reference $assembly > $out_dir/intermediate_deletions_HC.$date.vcf
exit_value=$?; if [ $exit_value != 0 ]; then echo "<< Error >>"; exit 1; else echo "<< Finished >>"; fi
rm $out_dir/IMSindel_overall_deletions.$date
rm $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD_RS_filtered.BP.intermediate
rm $out_dir/IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD_RS_filtered.BP.intermediate.temp
echo "<< Complete >>"
