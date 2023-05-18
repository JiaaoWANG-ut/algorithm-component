#!/bin/bash -l
#SBATCH -p small
#SBATCH -n 56
#SBATCH -N 1
#SBATCH -t 8:00:00
#SBATCH -J cp2k-HT
#SBATCH -V
#SBATCH -A CHE21006
#SBATCH -o ll_out

export OMP_NUM_THREAD=56

touch Started
/usr/local/bin/ibrun   /work2/08197/jiaao/frontera/apps/cp2k-7.1/exe/local/cp2k.popt -i *.inp
touch Finished
