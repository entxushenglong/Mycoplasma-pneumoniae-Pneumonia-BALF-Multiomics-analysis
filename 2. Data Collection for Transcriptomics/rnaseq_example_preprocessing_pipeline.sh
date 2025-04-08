{\rtf1\ansi\ansicpg936\cocoartf2821
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww19280\viewh15160\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #!/bin/bash\
\
# Step 1: Perform quality control using fastp\
fastp -i /path/to/rawdata/example/example_1.fq.gz -I /path/to/rawdata/example/example_2.fq.gz \\\
  -o /path/to/QC/example_1.fq.gz -O /path/to/QC/example_2.fq.gz \\\
  -j /path/to/QC/example.json -h /path/to/QC/example.html\
echo `date`' done'\
\
# Step 2: Create test files to check for rRNA contamination using Bowtie2\
cd /path/to/QC/\
gzip -dc /path/to/QC/example_1.fq.gz | head -n 40000 > /path/to/QC/example_1.test.fq && \\\
gzip -dc /path/to/QC/example_2.fq.gz | head -n 40000 > /path/to/QC/example_2.test.fq && \\\
bowtie2 -x /path/to/reference/mm10/rRNA \\\
  -1 /path/to/QC/example_1.test.fq -2 /path/to/QC/example_2.test.fq \\\
  --end-to-end --sensitive -p 8 --phred33 --no-mixed -X 600 \\\
  -S /path/to/QC/example.test.rRNA.sam 2> /path/to/QC/example.test.rRNA.stat.txt && \\\
rm /path/to/QC/example_1.test.fq /path/to/QC/example_2.test.fq /path/to/QC/example.test.rRNA.sam && \\\
maprate=$(tail -n 1 /path/to/QC/example.test.rRNA.stat.txt | awk '\{print $1\}' | awk -F '%' '\{print $1\}') && \\\
maprate=$(echo $maprate/1 | bc) && \\\
\
# Step 3: If rRNA contamination is above threshold, remove rRNA reads\
if [ $maprate -ge 10 ]; then\
  bowtie2 -x /path/to/reference/mm10/rRNA \\\
    -1 /path/to/QC/example_1.fq.gz -2 /path/to/QC/example_2.fq.gz \\\
    --end-to-end --sensitive -p 8 --phred33 --no-mixed -X 600 \\\
    --un-conc-gz /path/to/QC/example.unmap.gz \\\
    -S /path/to/QC/example.rRNA.sam 2>/path/to/QC/example.rRNA.stat.txt && \\\
  rm /path/to/QC/example.rRNA.sam && \\\
  mv /path/to/QC/example.unmap.1.gz /path/to/QC/example_1.fq.gz && \\\
  mv /path/to/QC/example.unmap.2.gz /path/to/QC/example_2.fq.gz\
fi\
\
# Step 4: Align reads to the reference genome using HISAT2\
hisat2 -p 2 --dta-cufflinks -x /path/to/reference/mm10/hisat2/mm10 \\\
  -1 /path/to/QC/example_1.fq.gz -2 /path/to/QC/example_2.fq.gz --no-unal -t | \\\
  samtools view -Sb > /path/to/Mapping/example.ori.bam\
\
# Step 5: Sort and index the BAM file\
samtools sort /path/to/Mapping/example.ori.bam -O bam -o /path/to/Mapping/example.sort.bam && \\\
samtools index /path/to/Mapping/example.sort.bam && \\\
echo `date`' done'\
\
# Step 6: Count reads mapped to genes using HTSeq\
python -m HTSeq.scripts.count -m union -s no -t exon -f bam \\\
  /path/to/Mapping/example.sort.bam \\\
  /path/to/Quantification/all_transcripts.gtf > \\\
  /path/to/Quantification/example.genes.readcount\
echo `date`' done'\
\
# Step 7: Calculate RPKM (Reads Per Kilobase of transcript, per Million mapped reads)\
perl /path/to/Pipeline/bin/calRPKM.pl \\\
  /path/to/Quantification/example.genes.readcount \\\
  /path/to/Quantification/all_transcripts.gtf \\\
  /path/to/Quantification/example.genes.fpkm}