#PBS -q batch
#PBS -N STAR_align_all_tis
#PBS -l nodes=1:ppn=30
#PBS -l mem=90gb
#PBS -l walltime=36:00:00
#PBS -M jpm73279@uga.edu
#PBS -m ae
#PBS -o /scratch/jpm73279/04.lncRNA/02.Analysis/31.Generate_STAR_Alignment 
#PBS -e /scratch/jpm73279/04.lncRNA/02.Analysis/31.Generate_STAR_Alignment 
#PBS -j oe


cd $PBS_O_WORKDIR

ml Bowtie2
ml Anaconda3
ml STAR

source activate snake_make
snakemake -rs STAR_align_RNA.snake --cores 30
