#! /bin/bash

dir=~/Data/2019-10-01-Fabian-script/salmonella
gtf=~/Data/2019-10-01-Fabian-script/salmonella/index/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_sl1344.ASM21085v2.45.gtf
### Creating required directories
mkdir $dir/cutadapt_output $dir/alignment  $dir/non_unique_count $dir/unique_count

### Assigning the sequence of primers
GAT27dT=GTGAGTGATGGTTGAGGATGTGTGGAGNNNNN
GAT27dTRC=NNNNNCTCCACACATCCTCAACCATCACTCAC
GAT216N3G=GATGGTTGAGGATGTGTGGAGNNNNNN
GAT216N3GRC=NNNNNNCTCCACACATCCTCAACCATC
GAT27PCR=GTGAGTGATGGTTGAGGATGTGTGGAG
GAT27PCRRC=CTCCACACATCCTCAACCATCACTCAC

### Running Cutadapt
for i in {1..235}
do

  ##### Here we retrive the name of the fastq files to use it as an input for cutadapt run
  read1="$(find $dir/source -name '*_R1.fq.gz' -printf "%f\n" | sort | sed -n "$i"p)"
  read2="$(find $dir/source -name '*_R2.fq.gz' -printf "%f\n" | sort | sed -n "$i"p)"
  report="$(find $dir/source -name '*_R1.fq.gz' -printf "%f\n" | sort | sed -n "$i"p | awk '{gsub("_R1.fq.gz", "", $0); print}')"
  ##### Running the cutadapt command
  cutadapt -j 40 -u 15 -U 15 -u -3 -U -3 -m 20 -q 20 -a "A{100}" -A "A{100}" \
    -g $GAT27dT -G $GAT27dT -g $GAT27dTRC -G $GAT27dTRC \
    -g $GAT216N3G -G $GAT216N3G -g $GAT216N3GRC -G $GAT216N3GRC \
    -g $GAT27PCR -G $GAT27PCR -g $GAT27PCRRC -G $GAT27PCRRC \
    -o $dir/cutadapt_output/$read1 -p $dir/cutadapt_output/$read2 $dir/source/$read1 $dir/source/$read2 \
    > $dir/cutadapt_output/$report.txt

done


### Running STAR
for i in {1..235}
do

  read1="$(find $dir/cutadapt_output -name '*_R1.fq.gz' -printf "%f\n" | sort | sed -n "$i"p)"
  read2="$(find $dir/cutadapt_output -name '*_R2.fq.gz' -printf "%f\n" | sort | sed -n "$i"p)"
  report="$(find $dir/cutadapt_output -name '*_R1.fq.gz' -printf "%f\n" | sort | sed -n "$i"p | awk '{gsub("_R1.fq.gz", "", $0); print}')"
  ##### Running STAR
  STAR --runThreadN 40 \
	  --genomeDir $dir/index/starIndex \
	  --readFilesIn $dir/cutadapt_output/$read1 $dir/cutadapt_output/$read2 \
	  --readFilesCommand zcat \
	  --outFileNamePrefix $dir/alignment/$report \
	  --outSAMtype BAM Unsorted SortedByCoordinate \
	  --outWigType wiggle \
	  --outWigStrand Unstranded \
	  --limitBAMsortRAM 1283746182

done 



### Running featurecount for non unique reads
find $dir/alignment -name '*.sortedByCoord.out.bam' | parallel featureCounts -p -t exon -g gene_id -a $gtf -M -O --fraction -o {.}.non_unique.txt {}
### Running featurecount for unique reads
find $dir/alignment -name '*.sortedByCoord.out.bam' | parallel featureCounts -p -t exon -g gene_id -a $gtf -o {.}.unique.txt {}
