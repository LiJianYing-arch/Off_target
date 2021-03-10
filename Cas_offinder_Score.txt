/public/home/jyli/CRISPR/Method/Prediction_Offtarget/Cas_Offinder
>G1_sgRNA1
ATGGTTGCATCCTGCCTAAAAGG
>G1_sgRNA2
TTATACCATCGGACTTTGCTAGG	
>G2_sgRNA1
TTGGCTGCTCACTGTGTTTATGG
>G2_sgRNA2
AACCCCCCACCACTCCCCGGTGG
>G3_sgRNA1
GTCTACCATCAGTGTTCCCACGG
>G3_sgRNA2
TCATAATCCTCGATGATTTGTGG




~/software/CRISPR_soft/Cas-OFFinder/cas-offinder-master/cas-offinder GhsgRNA.input.txt C GhsgRNA.output.txt
awk '{OFS="\t"}{if($4 !~ /AA$/)print}' GhsgRNA.output.txt >GhsgRNA.output_PAM.txt


grep 'ATGGTTGCATCCTGCCTAAA' GhsgRNA.output_PAM.txt |awk '{OFS="\t"}{if($5 ~/+/)print $1,$2,$3+1,$3+20,$4,$5,$6}'  >G1_sgRNA1_S.txt
grep 'ATGGTTGCATCCTGCCTAAA' GhsgRNA.output_PAM.txt |awk '{OFS="\t"}{if($5 ~/-/)print $1,$2,$3+4,$3+23,$4,$5,$6}'  >G1_sgRNA1_A.txt
grep 'TTATACCATCGGACTTTGCT' GhsgRNA.output_PAM.txt |awk '{OFS="\t"}{if($5 ~/+/)print $1,$2,$3+1,$3+20,$4,$5,$6}'  >G1_sgRNA2_S.txt
grep 'TTATACCATCGGACTTTGCT' GhsgRNA.output_PAM.txt |awk '{OFS="\t"}{if($5 ~/-/)print $1,$2,$3+4,$3+23,$4,$5,$6}'  >G1_sgRNA2_A.txt
grep 'TTGGCTGCTCACTGTGTTTA' GhsgRNA.output_PAM.txt |awk '{OFS="\t"}{if($5 ~/+/)print $1,$2,$3+1,$3+20,$4,$5,$6}'  >G2_sgRNA1_S.txt
grep 'TTGGCTGCTCACTGTGTTTA' GhsgRNA.output_PAM.txt |awk '{OFS="\t"}{if($5 ~/-/)print $1,$2,$3+4,$3+23,$4,$5,$6}'  >G2_sgRNA1_A.txt
grep 'AACCCCCCACCACTCCCCGG' GhsgRNA.output_PAM.txt |awk '{OFS="\t"}{if($5 ~/+/)print $1,$2,$3+1,$3+20,$4,$5,$6}'  >G2_sgRNA2_S.txt
grep 'AACCCCCCACCACTCCCCGG' GhsgRNA.output_PAM.txt |awk '{OFS="\t"}{if($5 ~/-/)print $1,$2,$3+4,$3+23,$4,$5,$6}'  >G2_sgRNA2_A.txt
grep 'GTCTACCATCAGTGTTCCCA' GhsgRNA.output_PAM.txt |awk '{OFS="\t"}{if($5 ~/+/)print $1,$2,$3+1,$3+20,$4,$5,$6}'  >G3_sgRNA1_S.txt
grep 'GTCTACCATCAGTGTTCCCA' GhsgRNA.output_PAM.txt |awk '{OFS="\t"}{if($5 ~/-/)print $1,$2,$3+4,$3+23,$4,$5,$6}'  >G3_sgRNA1_A.txt
grep 'TCATAATCCTCGATGATTTG' GhsgRNA.output_PAM.txt |awk '{OFS="\t"}{if($5 ~/+/)print $1,$2,$3+1,$3+20,$4,$5,$6}'  >G3_sgRNA2_S.txt
grep 'TCATAATCCTCGATGATTTG' GhsgRNA.output_PAM.txt |awk '{OFS="\t"}{if($5 ~/-/)print $1,$2,$3+4,$3+23,$4,$5,$6}'  >G3_sgRNA2_A.txt

cat G1_sgRNA1_S.txt G1_sgRNA1_A.txt >G1_sgRNA1.txt
cat G1_sgRNA2_S.txt G1_sgRNA2_A.txt >G1_sgRNA2.txt
cat G2_sgRNA1_S.txt G2_sgRNA1_A.txt >G2_sgRNA1.txt
cat G2_sgRNA2_S.txt G2_sgRNA2_A.txt >G2_sgRNA2.txt
cat G3_sgRNA1_S.txt G3_sgRNA1_A.txt >G3_sgRNA1.txt
cat G3_sgRNA2_S.txt G3_sgRNA2_A.txt >G3_sgRNA2.txt



cut -f5 G1_sgRNA1.txt |tr '[:lower:]' '[:upper:]' >AP2_sgRNA1.offtarget.txt
cut -f5 G1_sgRNA2.txt |tr '[:lower:]' '[:upper:]' >AP2_sgRNA2.offtarget.txt
cut -f5 G2_sgRNA1.txt |tr '[:lower:]' '[:upper:]' >MYB44_sgRNA1.offtarget.txt
cut -f5 G2_sgRNA2.txt |tr '[:lower:]' '[:upper:]' >MYB44_sgRNA2.offtarget.txt
cut -f5 G3_sgRNA1.txt |tr '[:lower:]' '[:upper:]' >ARC_sgRNA1.offtarget.txt
cut -f5 G3_sgRNA2.txt |tr '[:lower:]' '[:upper:]' >ARC_sgRNA2.offtarget.txt

echo `seq 1 2037` |sed 's/ /\n/g' |sed 's/^/off_cal(ATGGTTGCATCCTGCCTAAAAGG,/g' |cut -d ',' -f1 >AP2_sgRNA1.ontarget.txt
echo `seq 1 2210` |sed 's/ /\n/g' |sed 's/^/off_cal(TTATACCATCGGACTTTGCTAGG,/g' |cut -d ',' -f1 >AP2_sgRNA2.ontarget.txt
echo `seq 1 2883` |sed 's/ /\n/g' |sed 's/^/off_cal(TTGGCTGCTCACTGTGTTTATGG,/g' |cut -d ',' -f1 >MYB44_sgRNA1.ontarget.txt
echo `seq 1 288` |sed 's/ /\n/g' |sed 's/^/off_cal(AACCCCCCACCACTCCCCGGTGG,/g' |cut -d ',' -f1 >MYB44_sgRNA2.ontarget.txt
echo `seq 1 1270` |sed 's/ /\n/g' |sed 's/^/off_cal(GTCTACCATCAGTGTTCCCACGG,/g' |cut -d ',' -f1 >ARC_sgRNA1.ontarget.txt
echo `seq 1 3653` |sed 's/ /\n/g' |sed 's/^/off_cal(TCATAATCCTCGATGATTTGTGG,/g' |cut -d ',' -f1 >ARC_sgRNA2.ontarget.txt

head -n 2037 PAM.txt >AP2_sgRNA1.PAM.txt
head -n 2210 PAM.txt >AP2_sgRNA2.PAM.txt
head -n 2883 PAM.txt >MYB44_sgRNA1.PAM.txt
head -n 288 PAM.txt >MYB44_sgRNA2.PAM.txt
head -n 1270 PAM.txt >ARC_sgRNA1.PAM.txt
head -n 3653 PAM.txt >ARC_sgRNA2.PAM.txt

paste -d ',' AP2_sgRNA1.ontarget.txt AP2_sgRNA1.offtarget.txt AP2_sgRNA1.PAM.txt > AP2_sgRNA1Score_Sub.txt
paste -d ',' AP2_sgRNA2.ontarget.txt AP2_sgRNA2.offtarget.txt AP2_sgRNA2.PAM.txt > AP2_sgRNA2Score_Sub.txt
paste -d ',' MYB44_sgRNA1.ontarget.txt MYB44_sgRNA1.offtarget.txt MYB44_sgRNA1.PAM.txt > MYB44_sgRNA1Score_Sub.txt
paste -d ',' MYB44_sgRNA2.ontarget.txt MYB44_sgRNA2.offtarget.txt MYB44_sgRNA2.PAM.txt > MYB44_sgRNA2Score_Sub.txt
paste -d ',' ARC_sgRNA1.ontarget.txt ARC_sgRNA1.offtarget.txt ARC_sgRNA1.PAM.txt > ARC_sgRNA1Score_Sub.txt
paste -d ',' ARC_sgRNA2.ontarget.txt ARC_sgRNA2.offtarget.txt ARC_sgRNA2.PAM.txt > ARC_sgRNA2Score_Sub.txt

cat OT_sgRNAScore.pl AP2_sgRNA1Score_Sub.txt >AP2_sgRNA1.pl
cat OT_sgRNAScore.pl AP2_sgRNA2Score_Sub.txt >AP2_sgRNA2.pl
cat OT_sgRNAScore.pl MYB44_sgRNA1Score_Sub.txt >MYB44_sgRNA1.pl
cat OT_sgRNAScore.pl MYB44_sgRNA2Score_Sub.txt >MYB44_sgRNA2.pl
cat OT_sgRNAScore.pl ARC_sgRNA1Score_Sub.txt >ARC_sgRNA1.pl
cat OT_sgRNAScore.pl ARC_sgRNA2Score_Sub.txt >ARC_sgRNA2.pl

perl AP2_sgRNA1.pl |sort -k1,1 > AP2_sgRNA1.offtargetScore.txt
perl AP2_sgRNA2.pl |sort -k1,1 > AP2_sgRNA2.offtargetScore.txt
perl MYB44_sgRNA1.pl |sort -k1,1 > MYB44_sgRNA1.offtargetScore.txt
perl MYB44_sgRNA2.pl |sort -k1,1 > MYB44_sgRNA2.offtargetScore.txt
perl ARC_sgRNA1.pl |sort -k1,1 > ARC_sgRNA1.offtargetScore.txt
perl ARC_sgRNA2.pl |sort -k1,1 > ARC_sgRNA2.offtargetScore.txt

cut -f1-7 G1_sgRNA1.txt |tr '[:lower:]' '[:upper:]'|sort -k5,5 >AP2_sgRNA1_offtargetloci.txt
cut -f1-7 G1_sgRNA2.txt |tr '[:lower:]' '[:upper:]'|sort -k5,5 >AP2_sgRNA2_offtargetloci.txt
cut -f1-7 G2_sgRNA1.txt |tr '[:lower:]' '[:upper:]'|sort -k5,5 >MYB44_sgRNA1_offtargetloci.txt
cut -f1-7 G2_sgRNA2.txt |tr '[:lower:]' '[:upper:]'|sort -k5,5 >MYB44_sgRNA2_offtargetloci.txt
cut -f1-7 G3_sgRNA1.txt |tr '[:lower:]' '[:upper:]'|sort -k5,5 >ARC_sgRNA1_offtargetloci.txt
cut -f1-7 G3_sgRNA2.txt |tr '[:lower:]' '[:upper:]'|sort -k5,5 >ARC_sgRNA2_offtargetloci.txt

paste AP2_sgRNA1_offtargetloci.txt AP2_sgRNA1.offtargetScore.txt|sed 's/SCAFFOLD/scaffold/g' >AP2_sgRNA1.final.txt
paste AP2_sgRNA2_offtargetloci.txt AP2_sgRNA2.offtargetScore.txt|sed 's/SCAFFOLD/scaffold/g' >AP2_sgRNA2.final.txt
paste MYB44_sgRNA1_offtargetloci.txt MYB44_sgRNA1.offtargetScore.txt|sed 's/SCAFFOLD/scaffold/g' >MYB44_sgRNA1.final.txt
paste MYB44_sgRNA2_offtargetloci.txt MYB44_sgRNA2.offtargetScore.txt|sed 's/SCAFFOLD/scaffold/g' >MYB44_sgRNA2.final.txt
paste ARC_sgRNA1_offtargetloci.txt ARC_sgRNA1.offtargetScore.txt|sed 's/SCAFFOLD/scaffold/g' >ARC_sgRNA1.final.txt
paste ARC_sgRNA2_offtargetloci.txt ARC_sgRNA2.offtargetScore.txt|sed 's/SCAFFOLD/scaffold/g' >ARC_sgRNA2.final.txt

cut -f2-4 AP2_sgRNA1.final.txt >AP2_sgRNA1.final.tmp
cut -f2-4 AP2_sgRNA2.final.txt >AP2_sgRNA2.final.tmp
cut -f2-4 MYB44_sgRNA1.final.txt >MYB44_sgRNA1.final.tmp
cut -f2-4 MYB44_sgRNA2.final.txt >MYB44_sgRNA2.final.tmp
cut -f2-4 ARC_sgRNA1.final.txt >ARC_sgRNA1.final.tmp
cut -f2-4 ARC_sgRNA2.final.txt >ARC_sgRNA2.final.tmp

bedtools intersect  -a AP2_sgRNA1.final.tmp -b ~/Public_data/TM-1_Genome/TM-1.Chr.gene.loci -wb|sort -k8nr >AP2_sgRNA1.final.gene.txt
bedtools intersect  -a AP2_sgRNA2.final.tmp -b ~/Public_data/TM-1_Genome/TM-1.Chr.gene.loci -wb |sort -k8nr >AP2_sgRNA2.final.gene.txt
bedtools intersect  -a MYB44_sgRNA1.final.tmp -b ~/Public_data/TM-1_Genome/TM-1.Chr.gene.loci -wb |sort -k8nr >MYB44_sgRNA1.final.gene.txt
bedtools intersect  -a MYB44_sgRNA2.final.tmp -b ~/Public_data/TM-1_Genome/TM-1.Chr.gene.loci -wb|sort -k8nr  >MYB44_sgRNA2.final.gene.txt
bedtools intersect  -a ARC_sgRNA1.final.tmp -b ~/Public_data/TM-1_Genome/TM-1.Chr.gene.loci -wb|sort -k8nr  >ARC_sgRNA1.final.gene.txt
bedtools intersect  -a ARC_sgRNA2.final.tmp -b ~/Public_data/TM-1_Genome/TM-1.Chr.gene.loci -wb|sort -k8nr  >ARC_sgRNA2.final.gene.txt












