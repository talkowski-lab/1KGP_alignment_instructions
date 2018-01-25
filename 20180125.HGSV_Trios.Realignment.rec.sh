PCR-free Illumina high coverage (75X) sequences of three trio: CHS, PUR, YRI, that have been produced by the HGSVC have been down sampled to ~30X and realigned to hs37d5 

working directory on erisone: 
```
/data/talkowski/Samples/1000Genomes/HGSV_Illumina_Alignment_GRCh38
```

Step1: extract fastqs from HGSV alignment in cram format:
```
bsub -q medium -J HG00512 -o HG00512.log -sla miket_sc "samtools fastq -1 HG00512.1.fq -2 HG00512.2.fq HG00512.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram"
bsub -q medium -J HG00513 -o HG00513.log -sla miket_sc "samtools fastq -1 HG00513.1.fq -2 HG00513.2.fq HG00513.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram"
bsub -q medium -J HG00514 -o HG00514.log -sla miket_sc "samtools fastq -1 HG00514.1.fq -2 HG00514.2.fq HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram"
bsub -q medium -J HG00731 -o HG00731.log -sla miket_sc "samtools fastq -1 HG00731.1.fq -2 HG00731.2.fq HG00731.alt_bwamem_GRCh38DH.20150715.PUR.high_coverage.cram"
bsub -q medium -J HG00732 -o HG00732.log -sla miket_sc "samtools fastq -1 HG00732.1.fq -2 HG00732.2.fq HG00732.alt_bwamem_GRCh38DH.20150715.PUR.high_coverage.cram"
bsub -q medium -J HG00733 -o HG00733.log -sla miket_sc "samtools fastq -1 HG00733.1.fq -2 HG00733.2.fq HG00733.alt_bwamem_GRCh38DH.20150715.PUR.high_coverage.cram"
bsub -q medium -J NA19238 -o NA19238.log -sla miket_sc "samtools fastq -1 NA19238.1.fq -2 NA19238.2.fq NA19238.alt_bwamem_GRCh38DH.20150715.YRI.high_coverage.cram"
bsub -q medium -J NA19239 -o NA19239.log -sla miket_sc "samtools fastq -1 NA19239.1.fq -2 NA19239.2.fq NA19239.alt_bwamem_GRCh38DH.20150715.YRI.high_coverage.cram"
bsub -q medium -J NA19240 -o NA19240.log -sla miket_sc "samtools fastq -1 NA19240.1.fq -2 NA19240.2.fq NA19240.alt_bwamem_GRCh38DH.20150715.YRI.high_coverage.cram"
```

Step2: downsample fastq to 30X

 - downsample read 1 from 75X to 30X:
```
/data/talkowski/xuefang/local/src/seqtk/seqtk sample -s100  HG00512.1.fq .4 > HG00512.1.30X_S100.fq  
```
 - extract read names from read 1
```
grep '@' HG00512.1.30X_S100.fq  > HG00512.30X_S100.names
```
 - extract reads with the same name from read2
```
/data/talkowski/xuefang/local/src/seqtk/seqtk subseq HG00512.2.fq HG00512.30X_S100.names > HG00512.2.30X_S100.fq  
```


