Illumina sequences alignment to hs37d5

## download [reference genome](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/):

```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/*
```
### index reference sequences
```
samtools faidx hs37d5.fa
bwa index hs37d5.fa
```

## download other sources for alignment:
### 1. Indels for realignment
 - [High quality, experiment-validated indel set produced by Devine and Mills.](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf.gz)
```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf.gz
```

### 2. SNPs for recalibration
 - [dbSNP 135](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.dbsnp.build135.snps.sites.vcf.gz)
```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.dbsnp.build135.snps.sites.vcf.gz
```


### 3. Align with bwa mem
```
bwa mem  -t 1 -B 4 -O 6 -E 1 -M  ../reference/hs37d5.fa sample.1.fq sample.2.fq | samtools view -1 - > sample.hs37d5.bam
```

