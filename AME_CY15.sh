ls /home/yasue/bigdata/yasue/PCCOfPCC/Motif/MultiFasta/Target > /home/yasue/bigdata/yasue/PCCOfPCC/Motif/MultiFasta/allfasta.txt
for i in `cat /home/yasue/bigdata/yasue/PCCOfPCC/Motif/MultiFasta/allfasta.txt`
do
ame --o /home/yasue/bigdata/yasue/PCCOfPCC/Motif/MultiFasta/Results/${i}_results --control /home/yasue/bigdata/yasue/PCCOfPCC/Motif/MultiFasta/Control/control_${i} /home/yasue/bigdata/yasue/PCCOfPCC/Motif/MultiFasta/Target/${i} /home/yasue/bigdata/yasue/MEME/Databases/motif_databases/CIS-BP/Arabidopsis_thaliana.meme
done
