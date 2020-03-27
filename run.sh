#! /bin/bash

dir=~/evafadar/2020-02-21-Fabian_revision_order_20-0034
gtf=~/evafadar/2020-02-21-Fabian_revision_order_20-0034/index/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_sl1344.ASM21085v2.45.gtf

### Creating required directories
mkdir $dir/cutadapt_output $dir/alignment  $dir/non_unique_count $dir/unique_count

### Assigning the sequence of primers
GAT27dT=GTGAGTGATGGTTGAGGATGTGTGGAGNNNNN
GAT27dTRC=NNNNNCTCCACACATCCTCAACCATCACTCAC
GAT216N3G=GATGGTTGAGGATGTGTGGAGNNNNNN
GAT216N3GRC=NNNNNNCTCCACACATCCTCAACCATC
GAT27PCR=GTGAGTGATGGTTGAGGATGTGTGGAG
GAT27PCRRC=CTCCACACATCCTCAACCATCACTCAC
Nextra=CTGTCTCTTATA
NextraRC=TATAAGAGACAG

### Running Cutadapt
for i in {1..131}
do

##### Here we retrive the name of the fastq files to use it as an input for cutadapt run
read1="$(find /vol/data/gmak/runs/analyzed/novaseq/200220_A00278_0129_AHFKCGDRXX/HIRI_SIGA/ -name '*_R1_001.fastq.gz' -printf "%f\n" | sort | sed -n "$i"p)"
read2="$(find /vol/data/gmak/runs/analyzed/novaseq/200220_A00278_0129_AHFKCGDRXX/HIRI_SIGA/ -name '*_R2_001.fastq.gz' -printf "%f\n" | sort | sed -n "$i"p)"
report="$(find /vol/data/gmak/runs/analyzed/novaseq/200220_A00278_0129_AHFKCGDRXX/HIRI_SIGA/ -name '*_R1_001.fastq.gz' -printf "%f\n" | sort | sed -n "$i"p | awk '{gsub("_R1_001.fastq.gz", "", $0); print}')"
##### Running the cutadapt command
cutadapt -j 20 -u 20 -U 20 -u -3 -U -3 -m 20 -q 30 -a "A{100}" -A "A{100}" \
    -g $GAT27dT -G $GAT27dT -g $GAT27dTRC -G $GAT27dTRC \
    -g $GAT216N3G -G $GAT216N3G -g $GAT216N3GRC -G $GAT216N3GRC \
    -g $GAT27PCR -G $GAT27PCR -g $GAT27PCRRC -G $GAT27PCRRC \
    -a $Nextra -A $Nextra -a $NextraRC -A $NextraRC \
    -o $dir/cutadapt_output/$read1 -p $dir/cutadapt_output/$read2 /vol/data/gmak/runs/analyzed/novaseq/200220_A00278_0129_AHFKCGDRXX/HIRI_SIGA/$read1 /vol/data/gmak/runs/analyzed/novaseq/200220_A00278_0129_AHFKCGDRXX/HIRI_SIGA/$read2 \
    > $dir/cutadapt_output/$report.txt

done


### Running STAR
for i in {1..131}
do

read1="$(find $dir/cutadapt_output -name '*_R1_001.fastq.gz' -printf "%f\n" | sort | sed -n "$i"p)"
read2="$(find $dir/cutadapt_output -name '*_R2_001.fastq.gz' -printf "%f\n" | sort | sed -n "$i"p)"
report="$(find $dir/cutadapt_output -name '*_R1_001.fastq.gz' -printf "%f\n" | sort | sed -n "$i"p | awk '{gsub("_R1_001.fastq.gz", "", $0); print}')"
##### Running STAR
STAR --runThreadN 20 \
     --genomeDir $dir/index/starIndex \
     --readFilesIn $dir/cutadapt_output/$read1 $dir/cutadapt_output/$read2 \
     --readFilesCommand zcat \
     --outFileNamePrefix $dir/alignment/$report \
     --outSAMtype BAM Unsorted SortedByCoordinate \
     --outWigType wiggle \
     --outWigStrand Unstranded \

done 



### Running featurecount for non unique reads
find $dir/alignment -name '*ed.out.bam' | parallel featureCounts -p -t exon -g gene_id -a $gtf -M -O --fraction -o {.}.non_unique.txt {}
### Running featurecount for unique reads
find $dir/alignment -name '*ed.out.bam' | parallel featureCounts -p -t exon -g gene_id -a $gtf -o {.}.unique.txt {}
