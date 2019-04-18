## **_Joint-call recovery & filtering of intermediate-sized deletions_**

A workflow that generates a list of high-confidence intermediate-sized deletions following indel calling using IMSindel software.

### **Required:**

IMSindel (https://github.com/NCGG-MGC/IMSindel)  
Python (ver. 3.4 or higher)  
Pandas module for Python (ver. 0.23.4 or higher)  
R (ver. 3.5.1 or higher)  

---

### **Input files:**

*.out files from IMSindel run on each individual sample.

---

### **Notes before starting:**

Annotation files used are mainly taken from, or modified from the UCSC genome browser database (https://genome.ucsc.edu/index.html). Please prepare the following files and place them into the "required_files" directory before running the scripts. Ensure that you **_rename the files_** as shown. 

**1. hg19_genes_ucsc.txt**

- This file can be downloaded from the UCSC genome browser site (https://genome.ucsc.edu/index.html) under the Table Browser utility. Ensure group="Genes and Gene Predictions"; track="NCBI RefSeq"; table="UCSC RefSeq (refGene)"; assembly="Feb.2009 (GRCh37/hg19)".  


**2. hg19_simple_repeats_ucsc.txt**

- This file can be downloaded from the UCSC genome browser site (https://genome.ucsc.edu/index.html) under the Table Browser utility. Ensure group="Repeats"; track="Simple Repeats"; assembly="Feb.2009 (GRCh37/hg19)".  


**3. hg19_repeatmasker_ucsc.txt**

- This file can be downloaded from the UCSC genome browser site (https://genome.ucsc.edu/index.html) under the Table Browser utility. Ensure group="Repeats"; track="RepeatMasker"; assembly="Feb.2009 (GRCh37/hg19)".  


**4. hg19_telomere_centromere.txt**

- This file is provided. The file is modified after the original data was downloaded from the UCSC genome browser site (https://genome.ucsc.edu/index.html) under the Table Browser utility, with group="All Tables"; database=hg19; table="gap".  


**5. hg19_repeatmasker_ucsc_BED.txt**

- This file is provided as a .zip file. The file is modified from the "hg19_repeatmasker_ucsc.txt" file. Please expand the file before use.  


**6. bam_files_locations.txt**

- Prepare a tab-delimited text file containing sample names and locations of the corresponding .bam files of the samples.
First column should contain the sample name (e.g. "Sample1"), followed by the location of the .bam file in the second column (e.g. "/path/to/directory/containing/Sample1.bam").  

- Ensure that both .bam and .bai index files are located in the same directory/location. Name the file **_"bam_files_locations.txt"_** and place into _"required_files"_ directory.  

**Note: The "bam_files_locations.txt" file provided is only a placeholder file as an example of the file format.  
The .bam files listed (together with the .bai files) can be obtained from the 1000Genomes Projects ftp server from the URLs listed in the example file.**  
Please provide a file that **_corresponds to your own data_** should you wish to process your own samples!


**Example files:**

- Example files for a number of samples from the 1000Genome Project are provided in the directory "Example_files".  
- These files were produced from the outputs of indel calling using IMSindel software and processing the individual output files (as per the script <process_IMSindel_results.sh>) and extracting information for chromosome 21.   
- Please use these files if you would like to test the <generate_high_confidence_deletions.sh> script.  

---

### **Usage:**
1. Make an overall deletions list from results of IMSindel indel call results for each sample run.

- Use the script <process_IMSindel_results.sh> to process the .out files of IMSindel results.  
- Ensure that each samples' .out files are contained in a single individual directory specific to the sample. E.g. _/path/to/results/sample1_results_ and _/path/to/results/sample2_results_ ...etc.  
- It is recommended to have all samples' deletions lists within a single directory, i.e. /path/to/all/samples_deletions/  
- Deletions from the IMSindel call for each sample will be extracted from the results of the IMSindel run and collected into a single deletions list per sample. 
```
$ sh process_IMSindel_results.sh <path-to-individual-sample-IMSindel-results> <samplename(e.g. NA18943)> <path-to-output-directory-for-all-samples-deletion-lists>
```

2. Generate high-confidence deletions list from deletions of all the samples.

- Deletions from individual samples merged and high-confidence deletion list generated by filtering and joint-call recovery steps.  
- Use script <generate_high_confidence_deletions.sh> to generate list of high-confidence deletions and the VCF file.  
```
$ sh generate_high_confidence_deletions.sh <directory-containing-all-samples-deletion-lists> <output-results-directory>
```
---

### **Output files:**

The following file will be output for each sample for the <process_IMSindel_results.sh> script:

i)    **samplename_overall.results.chrm.deletions.date :** Overall deletions file for each sample

- This file contains the deletions extracted from the results of running IMSindel for each sample. Columns are tab-separated values.  
- "sample" column consists of the sample name as specified when making the overall merged file.  
- "contig", "sttpos" and "endpos" columns shows the contig in which the deletion was called from during IMSindel calling, as well as the start and end breakpoints.  
- "chr", "chrm_stt" and "chrm_edd" columns show the chromosome, start and end breakpoints, respectively, for the deletion. 
- All other columns are outputs from IMSindel.  

The following files will be output from running the <generate_high_confidence_deletions.sh> script:

i)     **IMSindel.overall.deletions.merged.date :** Merged deletions list from multiple samples  
  
ii)    **IMSindel.deletions.annotated.date :** Annotated deletions list 
  
iii)   **IMSindel.deletions.annotated.filtered.date :** Deletions list filtered by annotations 
  
iv)    **IMSindel.deletions.annotated.filtered.joint_recovery.date :** Deletions list for joint-call recovery use 
  
v)     **IMSindel.deletions.annotated.joint_recovery.filtered.date :** Deletions list after joint-call recovery 
  
vi)    **IMSindel.deletions.annotated.joint_recovery.HWE_filtered.date :** Deletions list after filtering by HWE 
  
vii)   **IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD_RS.date :** Deletions list for use to filter by read   depth of deletion region and average deletion breakpoints quality scores 
  
viii)  **IMSindel.deletions.annotated.joint_recovery.HWE_filtered.RD_RS_filtered.date :** Deletions list after filtering with  deletion read depth and average breakpoints quality score and with deletions flanking regions annotated with repeat-masked regions 
  
ix)    **intermediate_deletions_HC.date :** High-confidence deletions list after selecting for candidates with flanking regions not overlapped by same-class repeats 
  
x)     **intermediate_deletions_HC.date.vcf :** High-confidence deletions in VCF format 

The columns in the outputs files are tab-separated.

- "chr" , "stt" and "edd" columns indicate the chromosome, start breakpoint, and end breakpoint, respectively, for the detected deletions. 

- "gene_name" shows the gene/gene transcript in which the deletion is located within or overlaps. If the deletion lies completely within the gene, it is annotated as "within", while the "partial" annotation indicates the deletion overlaps a gene incompletely.  

- "gene_annotation" indicates whether a deletion overlaps a gene exon and the proprtion of overlap.  

- "centromere/telomere" column indicated if a deletion is located in a centromeric or telomeric region of the chromosome.

- "simple_repeats" and "repeatmasker" columns indicate the simple repeats and repeat-masked regions, respectively, in which the deletion may locate within or overlap, as well as the proportion of overlap.  

- "forward_flanking_region" and "reverse_flanking_region" columns indicate if the front and reverse flanking 100bp region of the deletion overlaps any simple repeat or repeat-masked region.  

- "Samples_num" column indicate the number of samples in which the deletion was detected in. 

- "Deletion_classes" show the numbers of types that the deletions were detected as in the samples. "SID" indicate the number of samples in which the deletion was called as "within-read" type. "B" and "F" indicate the number of samples in which the deletion as called using forward reads only, and reverse reads only, respectively. "LD" indicate the number of samples in which the deletion was called using both forward and reverse reads. "LD_sr" indicate the number of samples in which there were at least 5 forward and reverse reads that support the deletion call was observed. "LI" indicate the number of samples in which the deletion was called by IMSindel as "long insertion" and is used for filtering.  

- "obs_Homo", "obs_hetero", "obs_HomoWT" columns indicate the number of samples in which the deletion was observed to be homozygous, heterozygous or homozygous-reference allele, respectively. The columns <exp_Homo>, <exp_Hetero>, and <exp_HomoWT> columns indicate the expected numbers of homozygous, heterozygous, and homozygous-reference allele deletions. These values are used for calculating HWE p-values.  

- "P_HWE" column shows the HWE p-value calculated for the deletion.  

- "Avg_forward_region_readdepth", "Avg_deletion_readdepth", and "Avg_reverse_region_readdepth" columns show the average read depths calculated for the deletion region and its' flanking regions.  

- "Avg_forward_reads_softclip_score" and "Avg_reverse_reads_softclip_score" columns indicate the average quality scores of softclipped bases at the deletions' forward and reverse breakpoints, respectively.  

- "del_front_rep" and "del_rear_rep" columns indicate the repeat-masked regions annotation of the deletions' breakpoints for use in determining high-confidence candidates.  

- "samplenameX" columns contain information about the samples that the deletions were detected in. If the deletion was not detected in a particular sample, the information field is left as "-".  

The VCF file of the high-confidence deletions follows VCFv4.1 format. Deletion allele types are indicated as "0/0" for non-deletions (reference allele), "0/1" for heterozygous deletions, and "1/1" for homozygous deletions. These allele statuses are taken from the results of IMSindel indel calling.  

---

### **Parameter settings in configuration file:**

If you would like to change the parameters, please make changes to the provided config file. 

- merge_deletion_breakpoints_start: Size range (in basepairs) that deletions start breakpoints must be within to be considered as overlapping and collasped into a singe call. Default value is 10.

- merge_deletion_breakpoints_end: Size range (in basepairs) that deletions end breakpoints must be within to be considered as overlapping and collasped into a single call. Default value is 10.

- annotation_flank: Size (in basepairs) of the region adjacent to deletions breakpoints to be annotated with repeats. Used for filtering deletions. Default value is 100.

- simplerepeat_proportion_threshold: Threshold proportion of overlap between repeats (both simple repeats and repeat-masked regions) and deletion regions. Used for filtering deletions. Default value is 0.5.

- flanking_repeat_coverage_proportion_threshold: Threshold proportion of overlap between repeats (both simple repeats and repeat-masked regions) and the forward and reverse flanking regions of deletion regions. Used for filtering deletions. Default value is 0.5.

- LD_dels_support_reads_threshold: Threshold number of both forward- and reverse-type reads that support deletion call in a sample. Used for joint call recovery step. Default number is 5.

- SID_dels_number_threshold: Threshold number of samples which had deletion candidate of length less than 50bp called as "within-type". Used for joint call recovery step. Default number is 2.

- LD_dels_number_threshold: Threshold number of samples in which deletion was called using both forward- and reverse-type reads. Used for joint call recovery step. Default number is 2.

- LD_sr_number_threshold: Threshold number of samples in which deletion was seen to to be called by both forward- and reverse-type reads of at least the threshold number (LD_dels_support_reads_threshold). Used for joint call recovery step. Default number is 2.

- del_len_recovery_threshold: Size threshold of deletion (in basepairs) for use in categorizing deletion groups during joint call recovery step. Default value is 50.

- del_len_filter_threshold: Threshold of deletion size (in basepairs) for inclusion as intermediate-sized deletion. Default value is 30.

- HWE_exclusion_threshold: Hardy-Weinberg Equilibiurm p-value threshold for exclusion of deletion candidate. Default value is 0.0001.

- bamfile_extract_range: Size range (in basepairs) from the start and end breakpoints to extract from BAM file. Used for calculating the read-depth and average quality scores around breakpoints of deletion candidates. Default value is 1000.

- breakpoint_score_extract_range: Size range (in basepairs) around the start and end breakpoints of deletion candidates to use when calculating average quality scores. Default value is 110.

- readdepth_filter_threshold: Threshold value of read-depth for exclusion of deletion candidate. Default value is 150.

- readscore_filter_threshold: Threshold value of average quality score of the region around the deletion breakpoints. Default value is 15.

- reference: Source of reference used for the sequencing data. Value should be a string describing the reference source (e.g. NCBI37).  

- assembly: Source of genome assembly used for the sequencing data. Value should be a string describing the assembly (e.g.ncbi_build37.fa).

---

### **Performace:**

### **License:**

GPL

### **Contact:**
James Wong jing Hao (jh-wong@ddm.med.kyoto-u.ac.jp)
