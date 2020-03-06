for ID in s1 s3 s4 s15 s20 s23 s4 s62 s66 s75 s76 s77 s79 s84 s85 s91 s94
do
/usr/bin/java -Xmx30g -jar /public/home/cotton/software/GATK-3.1.1/GenomeAnalysisTK.jar -R /public/home/jyli/Public_data/TM-1_Genome/TM-1_genomeV01_2k.fa -T SelectVariants --variant /public/home/jyli/CRISPR/Off_Target_Project/VCF_file/gatk/${ID}.gatk.snp.vcf --concordance /public/home/jyli/CRISPR/Off_Target_Project/VCF_file/samtools/${ID}.samtools.snp.vcf -o ${ID}.gatk.samtools.commsnp.vcf 

#filter SNP Indel
MEANQUAL=`awk '{ if ($1 !~ /#/) { total += $6; count++ } }  END { print total/count }' ./${ID}.gatk.samtools.commsnp.vcf`

/usr/bin/java -Xmx5g -jar /public/home/cotton/software/GATK-3.1.1/GenomeAnalysisTK.jar  -R /public/home/jyli/Public_data/TM-1_Genome/TM-1_genomeV01_2k.fa -T VariantFiltration  --filterExpression "QD < 20.0 || ReadPosRankSum < -8.0 ||  FS > 10.0 || QUAL < $MEANQUAL "  --filterName LowQualFilter --variant ./${ID}.gatk.samtools.commsnp.vcf  --missingValuesInExpressionsShouldEvaluateAsFailing --logging_level ERROR -o ${ID}.concordance.flt1.snp.vcf

grep -v "Filter" ${ID}.concordance.flt1.snp.vcf >  ${ID}.concordance.filter1.snp.vcf
done


