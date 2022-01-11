for fq in /staging/nrr/ratib/T4/*.fq.gz
do
	fq=${fq##*/}
	sbatch --export=target_file=$fq breseq_T1clones.slurm
done