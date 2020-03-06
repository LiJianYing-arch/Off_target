#~/software/CRISPR_soft/Cas-OFFinder/cas-offinder-master/cas-offinder GhsgRNA.txt C GhsgRNA.out.txt
grep 'TTTAGGCAGGATGCAACCAT' GhsgRNA.output.txt |awk '{OFS="\t"}{if($5 ~/+/)print $1,$2,$3+1,$3+20,$4,$5,$6}'  >G1_sgRNA1_S.txt
grep 'TTTAGGCAGGATGCAACCAT' GhsgRNA.output.txt |awk '{OFS="\t"}{if($5 ~/-/)print $1,$2,$3+4,$3+23,$4,$5,$6}'  >G1_sgRNA1_A.txt
grep 'TTATACCATCGGACTTTGCT' GhsgRNA.output.txt |awk '{OFS="\t"}{if($5 ~/+/)print $1,$2,$3+1,$3+20,$4,$5,$6}'  >G1_sgRNA2_S.txt
grep 'TTATACCATCGGACTTTGCT' GhsgRNA.output.txt |awk '{OFS="\t"}{if($5 ~/-/)print $1,$2,$3+4,$3+23,$4,$5,$6}'  >G1_sgRNA2_A.txt
