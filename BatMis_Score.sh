cd /public/home/jyli/CRISPR/Method/Prediction_Offtarget/BatMis/NAG
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


grep 'G1_sgRNA1' sgRNA_PAMNAG_M5.txt |sort -k3,3 >AP2_sgRNA1_offtargetloci.txt
grep 'G1_sgRNA2' sgRNA_PAMNAG_M5.txt |sort -k3,3 >AP2_sgRNA2_offtargetloci.txt
grep 'G2_sgRNA1' sgRNA_PAMNAG_M5.txt |sort -k3,3 >MYB44_sgRNA1_offtargetloci.txt
grep 'G2_sgRNA2' sgRNA_PAMNAG_M5.txt |sort -k3,3 >MYB44_sgRNA2_offtargetloci.txt
grep 'G3_sgRNA1' sgRNA_PAMNAG_M5.txt |sort -k3,3 >ARC_sgRNA1_offtargetloci.txt
grep 'G3_sgRNA2' sgRNA_PAMNAG_M5.txt |sort -k3,3 >ARC_sgRNA2_offtargetloci.txt


paste AP2_sgRNA1_offtargetloci.txt AP2_sgRNA1.offtargetScore.txt |cut -f1-5,7|grep ':+'|sed 's/:+/\t/g'|awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$6+19,"+",$7}' >AP2_sgRNA1.final.txt
paste AP2_sgRNA1_offtargetloci.txt AP2_sgRNA1.offtargetScore.txt |cut -f1-5,7|grep ':-'|sed 's/:-/\t/g'|awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6+3,$6+22,"-",$7}' |sort -k4n >>AP2_sgRNA1.final.txt
paste AP2_sgRNA2_offtargetloci.txt AP2_sgRNA2.offtargetScore.txt |cut -f1-5,7|grep ':+'|sed 's/:+/\t/g'|awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$6+19,"+",$7}' >AP2_sgRNA2.final.txt
paste AP2_sgRNA2_offtargetloci.txt AP2_sgRNA2.offtargetScore.txt |cut -f1-5,7|grep ':-'|sed 's/:-/\t/g'|awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6+3,$6+22,"-",$7}' |sort -k4n >>AP2_sgRNA2.final.txt

paste MYB44_sgRNA1_offtargetloci.txt MYB44_sgRNA1.offtargetScore.txt |cut -f1-5,7|grep ':+'|sed 's/:+/\t/g'|awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$6+19,"+",$7}' >MYB44_sgRNA1.final.txt
paste MYB44_sgRNA1_offtargetloci.txt MYB44_sgRNA1.offtargetScore.txt |cut -f1-5,7|grep ':-'|sed 's/:-/\t/g'|awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6+3,$6+22,"-",$7}' |sort -k4n >>MYB44_sgRNA1.final.txt
paste MYB44_sgRNA2_offtargetloci.txt MYB44_sgRNA2.offtargetScore.txt |cut -f1-5,7|grep ':+'|sed 's/:+/\t/g'|awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$6+19,"+",$7}' >MYB44_sgRNA2.final.txt
paste MYB44_sgRNA2_offtargetloci.txt MYB44_sgRNA2.offtargetScore.txt |cut -f1-5,7|grep ':-'|sed 's/:-/\t/g'|awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6+3,$6+22,"-",$7}' |sort -k4n >>MYB44_sgRNA2.final.txt

paste ARC_sgRNA1_offtargetloci.txt ARC_sgRNA1.offtargetScore.txt |cut -f1-5,7|grep ':+'|sed 's/:+/\t/g'|awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$6+19,"+",$7}' >ARC_sgRNA1.final.txt
paste ARC_sgRNA1_offtargetloci.txt ARC_sgRNA1.offtargetScore.txt |cut -f1-5,7|grep ':-'|sed 's/:-/\t/g'|awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6+3,$6+22,"-",$7}' |sort -k4n >>ARC_sgRNA1.final.txt
paste ARC_sgRNA2_offtargetloci.txt ARC_sgRNA2.offtargetScore.txt |cut -f1-5,7|grep ':+'|sed 's/:+/\t/g'|awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$6+19,"+",$7}' >ARC_sgRNA2.final.txt
paste ARC_sgRNA2_offtargetloci.txt ARC_sgRNA2.offtargetScore.txt |cut -f1-5,7|grep ':-'|sed 's/:-/\t/g'|awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6+3,$6+22,"-",$7}' |sort -k4n >>ARC_sgRNA2.final.txt
