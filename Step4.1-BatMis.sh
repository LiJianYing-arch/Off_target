#~/bin/bin/build_index ~/Public_data/TM-1_Genome/BatMisdb/tm_1_v01c01.fa
/public/home/jyli/bin/bin/batman -q sgRNA.fa -g /public/home/jyli/Public_data/TM-1_Genome/BatMisdb/tm_1_v01c01.fa -n4 -mall -o sgRNA.bin -l ./log > ./log

/public/home/jyli/bin/bin/batdecode -i sgRNA.bin -g /public/home/jyli/Public_data/TM-1_Genome/BatMisdb/tm_1_v01c01.fa -L ./log -o sgRNA.tmp.txt > ./log

samtools view -S sgRNA.tmp.txt >sgRNA_Gosspium.fa.txt
#perl off_targetsite_select.pl sgRNA_NGG_mismatch_5.txt >sgRNA_PAMNGG_M5.txt
