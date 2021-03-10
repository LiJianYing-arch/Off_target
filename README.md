#### Whole genome sequencing reveals rare off-target mutations and considerable inherent genetic or/and somaclonal variations in CRISPR-Cas9-edited cotton plants

##### Overview

This document contain six key supplementary step and data processing scripts to describe the variations calling, filtering variations, merged variations, genome-wide prediction off-target sites, and detecting CRISPR off-target mutations. We note that our analysis of off-target mutations workflow is accessible online at http://cotton.hzau.edu.cn/EN/uploads/data/off-target/. 

##### Citation

Li J†, Manghwar H†, Sun L†, Wang P, Wang G, Sheng H, Zhang J, Liu H, Qin L, Rui H, Li B, Lindsey K, Daniell H, Jin S* and Zhang X. (2019) Whole genome sequencing reveals rare off-target mutations and considerable inherent genetic or/and somaclonal variations in CRISPR-Cas9-edited cotton plants.  Plant Biotechnol J. 2019;  17:858-868.

##### Steps

###### 1. Data access and trimming data

The raw data were downloaded from NCBI Sequence Read Archive (SRA) BioProject ID: PRJNA380842. The raw reads were filtered with Trimmomatic software1 (Version 0.32, MINLEN:75).

```bash
#!/bin/bash
#PBS -N TrimData
#PBS -l nodes=1:ppn=5
#PBS -o Trim.out
#PBS -e Trim.err
#PBS -q batch
#PBS -l walltime=100:0:0
cd /public/home/jyli/RawData/Off_Target_Projet/SRA
for ID in AP2_s1 AP2_s3 AP2_s4 AP2_s15 AP2_s20 AP2_s23 MYB44_s62 MYB44_s75 MYB44_s76 MYB44_s77 Ne_s65 Ne_s66 Ne_s67 WT_s195 WT_s199 WT_s79
do
java -jar /public/home/cotton/software/Trimmomatic/trimmomatic-0.32.jar PE -threads 5 -phred33 -trimlog ${ID}_trim.log ${ID}_R1_001.fastq.gz ${ID}_R2_001.fastq.gz  ./${ID}_1.clean.fastq.gz ./${ID}_1.unpaired.fastq.gz ./${ID}_2.clean.fastq.gz ./${ID}_2.unpaired.fastq.gz ILLUMINACLIP:/public/home/cotton/software/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:75
done
```



###### 2. SNPs/indels calling

The SNPs and indels calling through Samtools and the Genome Analysis Toolkit (GATK) software and filtering steps. This scripts mainly involve in reads mapping TM-1 genome using BWA software, discarding duplicates reads, generation the sorted of *bam files, further the SNP calling. For the indels calling, the BAM files were realigned using GATK with (-T RealignerTargetCreator, IndelRealigner) to get accurate indels.

```bash
#!/bin/bash
#PBS -N s1
#PBS -l nodes=1:ppn=10
#PBS -l walltime=480:00:00
#PBS -e errlog.out
#PBS -q batch

echo "Start at"
date

cd $PBS_O_WORKDIR
# bwa mem mapping ;
read1=/public/home/jyli/RawData/Off_Target_Projet/clean_data/s1_H2MHNDMXX_L1_1.clean.fq.gz
read2=/public/home/jyli/RawData/Off_Target_Projet/clean_data/s1_H2MHNDMXX_L1_2.clean.fq.gz
name=s1
bwa mem -M -t 10 -R '@RG\tID:${name}\tSM:${name}_SM' /public/home/jyli/Public_data/HZAU_PacBio_GhV1/SNP_INDEX/Ghirsutum_genome.fasta ${read1} ${read2} >${name}.sam

# sort
java -Xmx20g -jar /public/home/cotton/software/picard-tools/SortSam.jar INPUT=${name}.sam OUTPUT=${name}_srt.bam SORT_ORDER=coordinate >${name}.sam_srt 2>err.${name}.sam_srt

# mark duplicates
java -Xmx20g -jar /public/home/cotton/software/picard-tools/MarkDuplicates.jar INPUT=${name}_srt.bam OUTPUT=${name}_srt_redup.bam METRICS_FILE=metrics.txt > ${name}_redup 2>err.${name}_redup

# index
java -Xmx20g -jar /public/home/cotton/software/picard-tools/BuildBamIndex.jar INPUT=${name}_srt_redup.bam >${name}_index 2>${name}_index

# target regions realignment
/usr/bin/java -Xmx20g -jar /public/home/cotton/software/GATK-3.1.1/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /public/home/jyli/Public_data/HZAU_PacBio_GhV1/SNP_INDEX/Ghirsutum_genome.fasta -I ${name}_srt_redup.bam -o ${name}_forIndelRealigner.intervals -nt 10 -allowPotentiallyMisencodedQuals

/usr/bin/java -Xmx20g -jar /public/home/cotton/software/GATK-3.1.1/GenomeAnalysisTK.jar -T IndelRealigner -R /public/home/jyli/Public_data/HZAU_PacBio_GhV1/SNP_INDEX/Ghirsutum_genome.fasta -I ${name}_srt_redup.bam -targetIntervals ${name}_forIndelRealigner.intervals -o ${name}_realigned.bam -allowPotentiallyMisencodedQuals

# gatk SNP calling
/usr/bin/java -Xmx30g -jar /public/home/cotton/software/GATK-3.1.1/GenomeAnalysisTK.jar -allowPotentiallyMisencodedQuals -R /public/home/jyli/Public_data/HZAU_PacBio_GhV1/SNP_INDEX/Ghirsutum_genome.fasta -T UnifiedGenotyper -glm snp -I ${name}_srt_redup.bam -o ${name}.gatk.snp.vcf -nt 10 -stand_call_conf 30.0 -stand_emit_conf 0

#samtools SNP calling
/public/home/cotton/software/samtools/samtools mpileup -DSugf /public/home/jyli/Public_data/HZAU_PacBio_GhV1/SNP_INDEX/Ghirsutum_genome.fasta ${name}_srt_redup.bam |/public/home/cotton/software/samtools/bcftools/bcftools view -Ncvg - >${name}.samtools.vcf

#indel calling
/usr/bin/java -Xmx30g -jar /public/home/cotton/software/GATK-3.1.1/GenomeAnalysisTK.jar -allowPotentiallyMisencodedQuals -R /public/home/jyli/Public_data/HZAU_PacBio_GhV1/SNP_INDEX/Ghirsutum_genome.fasta -T UnifiedGenotyper -glm INDEL -o ${name}.gatk.indel.vcf -nt 10 -rf BadCigar -I ${name}_realigned.bam -metrics ${name}.gatk.indel.metrics
rm -f *.sam *._srt.bam
echo "End at:"
date
```



###### 3. Common variations and Filtering variations 

The common variations detected by GATK and Samtools using GATK software merged the (-T Selectvariations) at least depth 20× for each site was retained for further analysis. Then, the SNPs/indels were filtered with GATK parameters (QD < 20.0 && ReadPosRankSum < -8.0 && FS > 10.0 && QUAL < $MEANQUAL). 

```bash
for ID in s1 s3 s4 s15 s20 s23 s4 s62 s66 s75 s76 s77 s79 s84 s85 s91 s94
do
/usr/bin/java -Xmx30g -jar /public/home/cotton/software/GATK-3.1.1/GenomeAnalysisTK.jar -R /public/home/jyli/Public_data/TM-1_Genome/TM-1_genomeV01_2k.fa -T SelectVariants --variant /public/home/jyli/CRISPR/Off_Target_Project/VCF_file/gatk/${ID}.gatk.snp.vcf --concordance /public/home/jyli/CRISPR/Off_Target_Project/VCF_file/samtools/${ID}.samtools.snp.vcf -o ${ID}.gatk.samtools.commsnp.vcf 

#filter SNP Indel
MEANQUAL=`awk '{ if ($1 !~ /#/) { total += $6; count++ } }  END { print total/count }' ./${ID}.gatk.samtools.commsnp.vcf`

/usr/bin/java -Xmx5g -jar /public/home/cotton/software/GATK-3.1.1/GenomeAnalysisTK.jar  -R /public/home/jyli/Public_data/TM-1_Genome/TM-1_genomeV01_2k.fa -T VariantFiltration  --filterExpression "QD < 20.0 || ReadPosRankSum < -8.0 ||  FS > 10.0 || QUAL < $MEANQUAL "  --filterName LowQualFilter --variant ./${ID}.gatk.samtools.commsnp.vcf  --missingValuesInExpressionsShouldEvaluateAsFailing --logging_level ERROR -o ${ID}.concordance.flt1.snp.vcf

grep -v "Filter" ${ID}.concordance.flt1.snp.vcf >  ${ID}.concordance.filter1.snp.vcf
done

##Merged variations from Cas9-edited and WT/Ne plants
The each vcf files were merged using using customized Perl scripts or the VCF files were merged from different groups with VCFtools. Finally, we selected the variations present in three WT plants, where the negative plants have some genotype but differ from the Cas9-edited plants.
perl getMergedSNP.pl Sample.id > Sample_merged.SNP.txt
```



###### 5. Genome-wide prediction off-target sites and their scores

The potential off-target sites were predicted (PAM = NGG, NAG, NGA) using BatMis and Cas-Offinder software. The most off-target sites were predicted based on CRISPR-P web tools, according to the scoring system of off-target sites by a high-throughput analysis in mammalian cells (Hsu, P.D. et al. NBT, 2013). Genome-wide prediction off-target site using BatMis software. Running this code get the potential off-target site (include NGG, NGA, and others). Of a total 4,413 potential off-target sites (allowing ≤ 5 mismatches within the 20-bp sgRNA and 3-bp PAM sequences, at least one tool detecting，including the 3,296 (NGG), 410 (PAM: NAG) and 707 (PAM: NGA)) were ptredicted by BatMis (2,114) and Cas_offinder (3,763) tools. 

```bash
#BatMis alignment using the sgRNAs

~/bin/bin/build_index ~/Public_data/TM-1_Genome/BatMisdb/tm_1_v01c01.fa
/public/home/jyli/bin/bin/batman -q sgRNA.fa -g /public/home/jyli/Public_data/TM-1_Genome/BatMisdb/tm_1_v01c01.fa -n4 -mall -o sgRNA.bin -l ./log > ./log
/public/home/jyli/bin/bin/batdecode -i sgRNA.bin -g /public/home/jyli/Public_data/TM-1_Genome/BatMisdb/tm_1_v01c01.fa -L ./log -o sgRNA.tmp.txt > ./log
samtools view -S sgRNA.tmp.txt >sgRNA_Gosspium.fa.txt

perl off_targetsite_select.pl sgRNA_NGG_mismatch_5.txt >sgRNA_PAMNGG_M5.txt
```

Filtering the off-target site with ≤ 4 mismatch from BatMis results using Perl script, this PAM divide into NGG, NAG, and NGA.

```bash
Runing this code:
#select NGG Mis <= 4
perl off_target.pl ../sgRNA_HZAU_Mis5_NGG.fa.txt >sgRNA_HZAU_Mis5_NGG.fa.txt
#select NAG Mis <= 3
perl off_target2.pl ../sgRNA_HZAU_Mis5_NAG.fa.txt >sgRNA_HZAU_Mis5_NAG.fa.txt
#select NGA Mis <=3
perl off_target2.pl ../sgRNA_HZAU_Mis5_NGA.fa.txt >sgRNA_HZAU_Mis5_NGA.fa.txt
```

To assess the potential off-target site score according to the scoring system of off-target sites in mammalian cells. Running this script get the score of specific off-target site score: Genome-wide prediction off-target site using Cas-OFFinder and CRISPR-offinder software.

```bash
perl OT_sgRNAScore.pl 
./cas-offinder GhsgRNA.txt C GhsgRNA.out.txt
perl CRISPR-offinder_1.2.pl -input Gh_target.fa -pamseq NGG -pamori 3 -pamlen 20 -gc_min 20 -gc_max 80 -mismatches 5 -strand b -gd ~/Public_data/TM-1_Genome/Genome/
```



###### 6. Detecting prediction off-target mutations

The final high-quality SNPs/indels were overlapped with the potential off-target site region using Bedtools packages. On- and Off-target mutations visualization via Integrative Genomics Viewer.

```bash
samtools view -h s195_realigned.bam "A06:1364270-1368857"  |samtools view -bS - >s195.bam
samtools sort s195.bam s195.sort
samtools index s195.sort.bam s195.sort.bam.bai
samtools view A06:1366270-1366857.s195.bam |awk '{OFS="\t"; print ">"$1"\n"$10}' - >Gh_A06G0136.s195.fastq
```