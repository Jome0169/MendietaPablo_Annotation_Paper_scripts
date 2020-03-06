#PBS -q batch
#PBS -N Iso-seq-Align
#PBS -l nodes=1:ppn=30
#PBS -l mem=90gb
#PBS -l walltime=36:00:00
#PBS -M jpm73279@uga.edu
#PBS -m ae
#PBS -o /scratch/jpm73279/04.lncRNA/05.Maize_v5/01.analysis/10.Iso-seq_analysis 
#PBS -e /scratch/jpm73279/04.lncRNA/05.Maize_v5/01.analysis/10.Iso-seq_analysis 
#PBS -j oe


cd $PBS_O_WORKDIR

ml Anaconda3


source activate snake_make
snakemake -rs iso_seq.snake --cores 30 --use-conda
