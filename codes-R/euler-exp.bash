#PBS -N RInfer
#PBS -l ncpus=8
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M matheus.saldanha@usp.br

TIMEFORMAT="%E";

# load modules
module load gcc/4.9.2;
module load R/3.5.0;

cd /home/mathjs/ElfProbTET/codes-R;
R -f inference.r;
