Whole genome sequencing reveals rare off-target mutations and considerable inherent genetic or/and somaclonal variations in CRISPR-Cas9-edited cotton plants

Overview
This document contain six key supplementary step and data processing scripts to describe the variations calling, filtering variations, merged variations, genome-wide prediction off-target sites, and detecting CRISPR off-target mutations.
We note that our analysis of off-target mutations workflow is accessible online at http://cotton.hzau.edu.cn/EN/uploads/data/off-target/. 
The test data is available through the following reference. 
https://onlinelibrary.wiley.com/doi/full/10.1111/pbi.13020
Citation
Li J†, Manghwar H†, Sun L†, Wang P, Wang G, Sheng H, Zhang J, Liu H, Qin L, Rui H, Li B, Lindsey K, Daniell H, Jin S* and Zhang X. (2019) Whole genome sequencing reveals rare off-target mutations and considerable inherent genetic or/and somaclonal variations in CRISPR-Cas9-edited cotton plants.  Plant Biotechnology Journal. 17, pp. 858–868.

Steps
1. Data access and trimming data
The raw data were downloaded from NCBI Sequence Read Archive (SRA) BioProject ID: PRJNA380842. The raw reads were filtered with Trimmomatic software1 (Version 0.32, MINLEN:75).
TrimData.pbs

2. SNPs/indels calling
The SNPs and indels calling through Samtools2 and the Genome Analysis Toolkit (GATK)3 software and filtering steps. This scripts mainly involve in reads mapping TM-1 genome4 using BWA5 software, discarding duplicates reads, generation the sorted of *bam files, further the SNP calling. For the indels calling, the BAM files were realigned using GATK with (-T RealignerTargetCreator, IndelRealigner) to get accurate indels.
Running this script for SNPs/indels calling: Step1-variation calling.pbs

3. Common variations and Filtering variations 
The common variations detected by GATK and Samtools using GATK software merged the (-T Selectvariations) at least depth 20× for each site was retained for further analysis. Then, the SNPs/indels were filtered with GATK parameters (QD < 20.0 && ReadPosRankSum < -8.0 && FS > 10.0 && QUAL < $MEANQUAL).
Step2-FiterlSNP.sh

4. Merged variations from Cas9-edited and WT/Ne plants
The each vcf files were merged using using customized Perl scripts or the VCF files were merged from different groups with VCFtools.

eg1: perl getMergedSNP.pl Sample.id > Sample_merged.SNP.txt
eg2: Step3-vcftools.pbs
Finally, we selected the variations present in three WT plants, where the negative plants have some genotype but differ from the Cas9-edited plants.
eg：(WT1==WT2==WT3==Ne1==Ne2==Ne3) & (WT1 != Cas9-edited)
awk '{OFS="\t"}{if($2 == $3 && $3 == $4 && $4 == $5 && $5 == $6 && $6 == $7 && $8 != $2)print $1}' offtarget.indel.allele.txt > ./AP2/s1.CRISPR.indel.id

5. Genome-wide prediction off-target sites and their scores
The potential off-target sites were predicted (PAM = NGG, NAG, NGA) using BatMis and Cas-Offinder6 software. The most off-target sites were predicted based on CRISPR-P7 web tools, according to the scoring system of off-target sites by a high-throughput analysis in mammalian cells8.

5.1 Method1: Genome-wide prediction off-target site using BatMis software
Runing this code get the potential off-target site (include NGG, NGA, and others)
Step4.1-BatMis.sh
Filtering the off-target site with ≤ 4 mismatch from BatMis results using Perl script, this PAM divide into NGG, NAG, and NGA.

Runing this code:
#select NGG Mis <= 4
perl off_target.pl ../sgRNA_HZAU_Mis5_NGG.fa.txt >sgRNA_HZAU_Mis5_NGG.fa.txt
#select NAG Mis <= 3
perl off_target2.pl ../sgRNA_HZAU_Mis5_NAG.fa.txt >sgRNA_HZAU_Mis5_NAG.fa.txt
#select NGA Mis <=3
perl off_target2.pl ../sgRNA_HZAU_Mis5_NGA.fa.txt >sgRNA_HZAU_Mis5_NGA.fa.txt

To assess the potential off-target site score according to the scoring system of off-target sites in mammalian cells.
Running this script get the score of specific off-target site score
perl OT_sgRNAScore.pl
BatMis_Score.sh

5.2 Method2: Genome-wide prediction off-target site using Cas-OFFinder software
~/software/CRISPR_soft/Cas-OFFinder/cas-offinder-master/cas-offinder GhsgRNA.txt C GhsgRNA.out.txt
Cas_offinder_Score.sh

5.3 Method3: Genome-wide prediction off-target site using CRISPR_offinder software
perl CRISPR-offinder_1.2.pl -input Gh_target.fa -pamseq NGG -pamori 3 -pamlen 20 -gc_min 20 -gc_max 80 -mismatches 5 -strand b -gd ~/Public_data/TM-1_Genome/Genome/

Of a total 4,413 potential off-target sites (allowing ≤ 5 mismatches within the 20-bp sgRNA and 3-bp PAM sequences, at least one tool detecting，including the 3,296 (NGG), 410
(PAM: NAG) and 707 (PAM: NGA)) were ptredicted by BatMis (2,114) and Cas_offinder (3,763) tools (Supplementary Dataset). 

6. Detecting prediction off-target mutations
The final high-quality SNPs/indels were overlapped with the potential off-target site region using Bedtools packages9. On- and Off-target mutations visualization via Integrative Genomics Viewer (IGV).
eg: get the specific loci for IGV 
samtools view -h s195_realigned.bam "A06:1364270-1368857"  |samtools view -bS - >s195.bam
samtools sort s195.bam s195.sort
samtools index s195.sort.bam s195.sort.bam.bai
samtools view A06:1366270-1366857.s195.bam |awk '{OFS="\t"; print ">"$1"\n"$10}' - >Gh_A06G0136.s195.fastq

References
1.	AM, B., M, L. & B, U. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics 30, 2114-2120 (2014).
2.	Li, H., Handsaker, B., Wysoker, A., Fennell, T. & Ruan, J. The Sequence Alignment-Map format and SAMtools. Bioinformatics 25, 2078-2079 (2009).
3.	Mckenna, A. et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res 20, 1297-1303 (2010).
4.	Zhang, T. et al. Sequencing of allotetraploid cotton (Gossypium hirsutum L. acc. TM-1) provides a resource for fiber improvement. Nat. Biotechnol. 33, 531-537 (2015).
5.	Li, H. & Durbin, R. Fast and accurate short read alignment with Burrows–Wheeler transform. Bioinformatics 25, 1754-1760 (2009).
6.	Bae, S., Park, J. & Kim, J.S. Cas-OFFinder: a fast and versatile algorithm that searches for potential off-target sites of Cas9 RNA-guided endonucleases. Bioinformatics 30, 1473-1475 (2014).
7.	Liu, H. et al. CRISPR-P 2.0: An Improved CRISPR-Cas9 Tool for Genome Editing in Plants. Mol Plant 10, 530-532 (2017).
8.  Hsu, P.D. et al. DNA targeting specificity of RNA-guided Cas9 nucleases. Nat. Biotechnol. 31, 827-832 (2013).
9. 	Quinlan, A.R. & Hall, I.M. BEDTools: a ﬂexible suite of utilities for comparing genomic features. Bioinformatics 26, 841–842 (2010).
