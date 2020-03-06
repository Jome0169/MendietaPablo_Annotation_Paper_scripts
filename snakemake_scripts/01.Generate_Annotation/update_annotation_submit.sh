#PBS -q batch
#PBS -N H3K36me3_Snake
#PBS -l nodes=1:ppn=30
#PBS -l mem=90gb
#PBS -l walltime=36:00:00
#PBS -M jpm73279@uga.edu
#PBS -m ae
#PBS -o /scratch/jpm73279/04.lncRNA/02.Analysis/24.regenerate_expression_peaks
#PBS -e /scratch/jpm73279/04.lncRNA/02.Analysis/24.regenerate_expression_peaks
#PBS -j oe


cd $PBS_O_WORKDIR

ml Trimmomatic
ml Anaconda3
ml STAR
ml Bowtie2


source activate snake_make
snakemake -rs Generate_peak_lists.snake --cores 30 --use-conda

