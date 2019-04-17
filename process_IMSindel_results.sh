#!/bin/bash

path=$0
if [ "$path" == "${path%/*}" ]; then src=.; else src=${path%/*}; fi
in_dir=$1
samplename=$2
date=`date "+%Y%m%d"`
out_dir=$3
overall_suffix="_overall.results"
chrm_suffix="_overall.results.chrm"
samplename_suffix="samplename"
deletions_file_suffix="deletions"
if [ ! -d $out_dir ]; then mkdir $out_dir; fi

echo "Merging IMSindel outputs and extracting deletions... "
cat $in_dir/*.out |python $src/bin/imsindel_results_clean.py /dev/stdin > $out_dir/$samplename$overall_suffix.$date
python $src/bin/attach_samplename.py $out_dir/$samplename$overall_suffix.$date $samplename > $out_dir/$samplename$chrm_suffix.$date.$samplename_suffix
#python $src/bin/attach_samplename.py $out_dir/$samplename$chrm_suffix.$date $samplename > $out_dir/$samplename$chrm_suffix.$date.$samplename_suffix
python $src/bin/extract_deletions_IMSindel_results.py $out_dir/$samplename$chrm_suffix.$date.$samplename_suffix > $out_dir/$samplename$chrm_suffix.$deletions_file_suffix.$date
rm $out_dir/$samplename$chrm_suffix.$date.$samplename_suffix
rm $out_dir/$samplename$overall_suffix.$date

echo "<< Complete >>"

