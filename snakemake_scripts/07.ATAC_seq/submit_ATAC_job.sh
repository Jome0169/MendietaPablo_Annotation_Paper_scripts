#PBS -S /bin/bash
#PBS -q batch
#PBS -N ATAC_snake
#PBS -l nodes=1:ppn=40
#PBS -l mem=100gb
#PBS -l walltime=15:00:00
#PBS -M jpm73279@uga.edu
#PBS -m abe
#PBS -o /scratch/jpm73279/04.lncRNA/02.Analysis/05.ATAC_align 
#PBS -e /scratch/jpm73279/04.lncRNA/02.Analysis/05.ATAC_align 
#PBS -j oe


cd $PBS_O_WORKDIR

ml picard 
ml MACS2 
ml deepTools 
ml Bowtie2


echo
echo "Job ID: $PBS_JOBID"
echo "Queue:  $PBS_QUEUE"
echo "Cores:  $PBS_NP"
echo "Nodes:  $(cat $PBS_NODEFILE | sort -u | tr '\n' ' ')"
echo "mpirun: $(which mpirun)"
echo


source activate snake_make
snakemake -rs  Align_ATAC_call_peaks.snake --cores 40 --use-conda
