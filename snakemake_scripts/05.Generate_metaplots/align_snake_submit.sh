#PBS -q batch
#PBS -N align_CHIP_seq
#PBS -l nodes=1:ppn=30
#PBS -l mem=90gb
#PBS -l walltime=15:00:00
#PBS -M jpm73279@uga.edu
#PBS -m ae
#PBS -o /scratch/jpm73279/04.lncRNA/02.Analysis/23.generate_all_metaplots
#PBS -e /scratch/jpm73279/04.lncRNA/02.Analysis/23.generate_all_metaplots
#PBS -j oe


cd $PBS_O_WORKDIR

ml Bowtie2
ml Anaconda3
ml deepTools
ml BEDTools

source activate snake_make
snakemake -rs align_chip_seq_reads.snake --cores 30
bash clean.sh
