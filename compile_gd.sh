for i in /Users/ratib/Dropbox/breseq/breseq_T4clones/*.fq.gz
do
	cp $i/output/output.gd /Users/ratib/Desktop/breseq_T3clones/gd/${i##*/}.gd
done