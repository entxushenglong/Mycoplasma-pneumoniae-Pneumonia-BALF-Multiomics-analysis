{\rtf1\ansi\ansicpg936\cocoartf2821
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #!/bin/bash\
\
# Activate kneaddata environment\
source activate kneaddata &&\
\
# Create necessary directories for storing intermediate and final results\
mkdir -p temp/unzip/ &&\
mkdir -p temp/addsuffix/ &&\
mkdir -p temp/kneaddataresults &&\
mkdir -p temp/humann &&\
mkdir -p result/metaphlan4 &&\
mkdir -p result/humann3 &&\
\
# Use tail and cut to extract sample IDs from metadata.tsv and parallelize file processing using rush\
tail -n+10 result/metadata79.tsv | cut -f1 | rush -k -j 4 \\\
"gunzip -c '/path/to/rawdata/\{\}_1.fq.gz' > 'temp/unzip/\{\}_1.fastq'; \\\
gunzip -c '/path/to/rawdata/\{\}_2.fq.gz' > 'temp/unzip/\{\}_2.fastq'; \\\
sed '1~4 s/ 1:/.1:/;1~4 s/$/\\/1/' temp/unzip/\{\}_1.fastq > temp/addsuffix/\{\}_1.fastq; \\\
sed '1~4 s/ 2:/.1:/;1~4 s/$/\\/2/' temp/unzip/\{\}_2.fastq > temp/addsuffix/\{\}_2.fastq; \\\
rm temp/unzip/\{\}_1.fastq; \\\
rm temp/unzip/\{\}_2.fastq; \\\
kneaddata -i1 temp/addsuffix/\{\}_1.fastq -i2 temp/addsuffix/\{\}_2.fastq \\\
-o temp/kneaddataresults --output-prefix \{\} \\\
-db /path/to/kneaddata/human_db \\\
--remove-intermediate-output -v -t 6; \\\
rm temp/addsuffix/\{\}_1.fastq; \\\
rm temp/addsuffix/\{\}_2.fastq; \\\
cat temp/kneaddataresults/\{\}_paired_1.fastq temp/kneaddataresults/\{\}_paired_2.fastq > temp/kneaddataresults/\{\}.fastq; \\\
rm temp/kneaddataresults/\{\}_paired_1.fastq; \\\
rm temp/kneaddataresults/\{\}_paired_2.fastq; \\\
rm temp/kneaddataresults/\{\}_unmatched_1.fastq; \\\
rm temp/kneaddataresults/\{\}_unmatched_2.fastq; \\\
rm temp/kneaddataresults/\{\}_hg37dec_v0.1_bowtie2* "\
\
# Activate HUMAnN3 environment\
source activate /path/to/humann3/env &&\
echo "Activated HUMAnN3 environment" &&\
\
# Run HUMAnN3 for microbiome functional analysis\
tail -n+2 result/metadata.tsv | cut -f1 | rush -k -j 8 \\\
"humann --input temp/kneaddataresults/\{\}.fastq --output temp/humann/ \\\
--threads 4 --metaphlan-options '--bowtie2db /path/to/metaphlan4_db --index mpa_vOct22_CHOCOPhlAnSGB_202212 --offline'; \\\
mv temp/humann/\{\}_humann_temp/\{\}_metaphlan_bugs_list.tsv temp/humann/ ; \\\
rm -r temp/humann/\{\}_humann_temp ; \\\
\
echo '\{\} HUMAnN3 processing finished'" &&\
echo "HUMAnN3 processing completed" &&\
\
# Merge Metaphlan result tables\
/path/to/script/merge_metaphlan_tables.py temp/humann/*_metaphlan_bugs_list.tsv | \\\
  sed 's/_metaphlan_bugs_list//g' | tail -n+2 | sed '1 s/clade_name/ID/' | sed '2i #metaphlan4'> result/metaphlan4/taxonomy.tsv &&\
\
# Convert the merged Metaphlan result to the stamp format\
/path/to/script/metaphlan_to_stamp.pl result/metaphlan4/taxonomy.tsv \\\
  |sort -r | uniq > result/metaphlan4/taxonomy.spf &&\
\
# Filter out 'unclassified' entries\
grep -v 'unclassified' result/metaphlan4/taxonomy.spf > result/metaphlan4/taxonomy2.spf &&\
\
# Merge HUMAnN result tables\
humann_join_tables --input temp/humann \\\
  --file_name pathabundance \\\
  --output result/humann3/pathabundance.tsv &&\
\
# Clean up file names in the pathabundance table\
sed -i 's/_Abundance//g' result/humann3/pathabundance.tsv &&\
\
# Normalize the abundance table (relative abundance)\
humann_renorm_table \\\
  --input result/humann3/pathabundance.tsv \\\
  --units relab \\\
  --output result/humann3/pathabundance_relab.tsv &&\
\
# Split the normalized table into stratified tables\
humann_split_stratified_table \\\
  --input result/humann3/pathabundance_relab.tsv \\\
  --output result/humann3/  &&\
\
# Convert HUMAnN results to KEGG annotations\
for i in `tail -n+2 result/metadata.tsv|cut -f1`;do\
  humann_regroup_table \\\
    -i temp/humann/$\{i\}_genefamilies.tsv \\\
    -g uniref90_ko \\\
    -o temp/humann/$\{i\}_ko.tsv\
done &&\
\
# Join the KO tables into a single file\
humann_join_tables \\\
  --input temp/humann/ \\\
  --file_name ko \\\
  --output result/humann3/ko.tsv &&\
\
# Clean up KO table file names\
sed -i '1s/_Abundance-RPKs//g' result/humann3/ko.tsv &&\
\
# Split the KEGG KO table into stratified tables\
humann_split_stratified_table \\\
  --input result/humann3/ko.tsv \\\
  --output result/humann3/  &&\
\
# Count the number of lines in the KO tables\
wc -l result/humann3/ko* &&\
\
# Summarize KEGG abundances\
/path/to/script/summarizeAbundance.py \\\
  -i result/humann3/ko_unstratified.tsv \\\
  -m /path/to/kegg/KO1-4.txt \\\
  -c 2,3,4 -s ',+,+,' -n raw \\\
  -o result/humann3/KEGG &&\
\
# Count the number of lines in the KEGG tables\
wc -l result/humann3/KEGG* &&\
\
# Final message indicating all processes are completed\
echo "All processes completed"}