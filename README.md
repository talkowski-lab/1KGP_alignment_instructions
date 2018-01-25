# illumina platinum pedigree realignment
This repository contains files and commands used to realign platinum genome to the GRCh38_BSM.fa

## download fastq file of platinum genome:
downloading address of platinum genome were indexed in illumina_platinum_ped.sequence.index

#### Programs used
The Illumina Platinum pedigree data was aligned to the reference genome using the following programs:

1. [BWA-MEM](https://github.com/lh3/bwa/blob/master/bwakit/README.md) 
[bwakit-0.7.13](https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.13_x64-linux.tar.bz2/download)
2. [Samtools-1.3](http://www.htslib.org/doc/samtools.html)(https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2)
3. [picard-tools-2.1.1](https://github.com/broadinstitute/picard/releases/download/2.1.1/picard-tools-2.1.1.zip)
[picard-2.1.1 source code](https://github.com/broadinstitute/picard/archive/2.1.1.zip)
4. [biobambam2](https://github.com/gt1/biobambam2/releases)
5. [GATK-3.3-0](https://github.com/broadgsa/gatk-protected/tree/3.3)
6. [Cramtools-3.0](https://github.com/enasequence/cramtools/tree/cram3)
7. [ant](http://supergsego.com/apache//ant/binaries/apache-ant-1.9.6-bin.zip)
8. [java1.8](http://download.oracle.com/otn-pub/java/jdk/8u73-b02/jdk-8u73-linux-x64.tar.gz)
```
setup java1.8 as default for ant:
pwd:/nfs/turbo/dcmb-brainsom/technical/application/picard-2.1.1
JAVA_HOME=/nfs/turbo/dcmb-brainsom/technical/application/jdk1.8.0_73/ ANT_HOME=/nfs/turbo/dcmb-brainsom/technical/application/apache-ant-1.9.6 "$ANT_HOME/bin/ant")
```

The Illumina Platinum pedigree data was aligned to the reference genome using the following reference datasets:

#### 1. Indels for realignment 
   - Phase 3 biallelic indels from Shapeit2 lifted to GRCh38 using the NCBI remapper ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz
   - High quality, experiment-validated indel set produced by Devine and Mills. The coordinates were lifted to GRCh38 by Alison Meynert from IGMM in Edinburgh using CrossMap and the UCSC chain files, filtered out ref == alt cases ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz

#### 2. SNPs for recalibration 
   - dbSNP 142 in GRCh38 coordinates ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz


## Command lines
1. Download and unzip fastq file for NA12878:
```
wget ERRXXX.fastq.gz
gunzip ERRXXX.fastq.gz
```

2. Alignment at run level
```
bwa mem  -t 1 -B 4 -O 6 -E 1 -M -R $rg_string $reference_fasta_file $fastq_file(1) $fastq_file(2) | samtools view -1 - > $bam_file
```

3. Local realignment around known indels by GATK
```
java -jar picard.jar CreateSequenceDictionary R=reference.fa O=reference.dict
java $jvm_args -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference_fasta -o $intervals_file -known $known_indels_file(s) 
java $jvm_args -jar GenomeAnalysisTK.jar -T IndelRealigner -R $reference_fasta -I $bam_file -o $realigned_bam_file -targetIntervals $intervals_file -known $known_indels_file(s) -LOD 0.4 -model KNOWNS_ONLY -compress 0 --disable_bam_indexing
```

4.  Recalibrate base quality scores using known SNPs by GATK
```
java $jvm_args -jar GenomeAnalysisTK.jar -T BaseRecalibrator -nt 1 -l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -R $reference_fasta -o $recal_data.table -I $bam_file -knownSites $known_snps_from_dbSNP142
java $jvm_args -jar GenomeAnalysisTK.jar -T PrintReads -l INFO -R $reference_fasta -o $recalibrated_bam -I $bam_file -BQSR $recal_data.table --disable_bam_indexing
```

5. Mark Duplicates using BioBamBam2
```
bammarkduplicates I=$input_bam O=$output_bam index=1 rmdup=0
```

6. Creating lossless CRAM files using Cramtools
```
java cramtools-3.0.jar cram --input-bam-file $input_bam --output-cram-file $output_cram --capture-all-tags --ignore-tags OQ:CQ:BQ --preserve-read-names --lossy-quality-score-spec *8 --reference-fasta-file $reference_fasta
```




## A quick example for NA12878
#### 1. Download and unzip fastq file for NA12878
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
gunzip ERR194147_1.fastq.gz
```

#### 2. Align NA12878 at run level
######## 2a. subsplit fastq file into smaller files for faster alignment
```
python /scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/pbs/subsplit_fastq.py --input /scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/fastq_2/NA12878_1.fastq --size 10000000
```
########2b. align sub fastq files
```
bwa mem  -t 1 -B 4 -O 6 -E 1 -M -R "@RG\tID:ERR194147\tSM:NA12878\tCN:ILLUMINA\tPL:ILLUMINA\tDS:ERP001960" /nfs/turbo/dcmb-brainsom/technical/reference/GRCh38_bsm_reference_genome/GRCh38_BSM.fa /scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/fastq_sub/NA12878_1.sub1.fastq /scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/fastq_sub/NA12878_2.sub1.fastq |samtools view -1 - > /scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/alignment/NA12878.sub1.bam
```
######## 2c. merge aligned sub bam files
```
samtools merge /scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/alignment/NA12878.1.bam /scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/alignment/NA12878.sub1.bam /scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/alignment/NA12878.sub2.bam /scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/alignment/NA12878.sub3.bam
```

#### 3. Local realignment around known indels by GATK
######## 3a. Use CreateSequenceDictionary.jar from Picard to create a .dict file from a fasta file
```
/nfs/turbo/dcmb-brainsom/technical/application/jdk1.8.0_73/bin/java -jar /nfs/turbo/dcmb-brainsom/technical/application/picard-2.1.1/dist/picard.jar CreateSequenceDictionary R=/nfs/turbo/dcmb-brainsom/technical/reference/GRCh38_bsm_reference_genome/GRCh38_BSM.fa O=/nfs/turbo/dcmb-brainsom/technical/reference/GRCh38_bsm_reference_genome/GRCh38_BSM.dict
```
######## 3b. Use RealignerTargetCreator from GATK to create a intervals_file
```
/nfs/turbo/dcmb-brainsom/technical/application/jdk1.8.0_73/bin/java -Xmx8g -Djava.io.tmpdir=temp -jar /nfs/turbo/dcmb-brainsom/technical/application/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /nfs/turbo/dcmb-brainsom/technical/reference/GRCh38_bsm_reference_genome/GRCh38_BSM.fa -o /nfs/turbo/dcmb-brainsom/technical/application/support_files/Platinum_BSM_intervals_file.picard -known /nfs/turbo/dcmb-brainsom/technical/application/support_files/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.bsm.modifed.vcf -known /nfs/turbo/dcmb-brainsom/technical/application/support_files/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf
```
######## 3c. Local realignment around known indels by GATK
```
/nfs/turbo/dcmb-brainsom/technical/application/jdk1.8.0_73/bin/java -Xmx8g -Djava.io.tmpdir=temp2.NA12878.sorted.bam -jar /nfs/turbo/dcmb-brainsom/technical/application/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T IndelRealigner -R /nfs/turbo/dcmb-brainsom/technical/reference/GRCh38_bsm_reference_genome/GRCh38_BSM.fa -I /scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/alignment/NA12878.sorted.bam -o /scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/alignment/NA12878.sorted.realign.bam -targetIntervals /nfs/turbo/dcmb-brainsom/technical/application/support_files/Platinum_BSM_intervals_file.picard -known /nfs/turbo/dcmb-brainsom/technical/application/support_files/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.bsm.modifed.vcf -known /nfs/turbo/dcmb-brainsom/technical/application/support_files/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf -LOD 0.4 -model KNOWNS_ONLY -compress 0 --disable_bam_indexing
```

#### 4.  Recalibrate base quality scores using known SNPs by GATK
```
/nfs/turbo/dcmb-brainsom/technical/application/jdk1.8.0_73/bin/java -Xmx8g -Djava.io.tmpdir=temp3 -jar /nfs/turbo/dcmb-brainsom/technical/application/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T BaseRecalibrator -nt 1 -l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -R /nfs/turbo/dcmb-brainsom/technical/reference/GRCh38_bsm_reference_genome/GRCh38_BSM.fa -o /scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/alignment/NA12878.sorted.recal_data.table -I /scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/alignment/NA12878.sorted.realign.bam -knownSites /nfs/turbo/dcmb-brainsom/technical/application/support_files/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf
/nfs/turbo/dcmb-brainsom/technical/application/jdk1.8.0_73/bin/java -Xmx8g -Djava.io.tmpdir=temp3 -jar /nfs/turbo/dcmb-brainsom/technical/application/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T PrintReads -l INFO -R /nfs/turbo/dcmb-brainsom/technical/reference/GRCh38_bsm_reference_genome/GRCh38_BSM.fa -o /scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/alignment/NA12878.sorted.realign.recalibrated.bam -I /scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/alignment/NA12878.sorted.realign.bam -BQSR /scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/alignment/NA12878.sorted.recal_data.table --disable_bam_indexing
```

#### 5. Mark Duplicates using BioBamBam2
```
/nfs/turbo/dcmb-brainsom/technical/application/biobambam2-2.0.35-release-20160330111451-x86_64-etch-linux-gnu/bin/bammarkduplicates I=/scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/alignment/NA12878.sorted.realign.recalibrated.bam O=/mnt/EXT/Mills-scratch2/Xuefang/bsmn/platinum.genome.realignment/alignment/NA12878.sorted.realign.recalibrated.markdup.bam index=1 rmdup=0
```

#### 6. Check quality of finalized bam files
```
samtools flagstat /scratch/remills_flux/xuefzhao/bsmn/platinum.genome.realignment/alignment/NA12878.sorted.realign.recalibrated.markdup.bam
```
