#!/bin/bash
#SBATCH -J test
#SBATCH --comment=WRF
#SBATCH -N 2
#SBATCH --ntasks-per-node=32
#SBATCH -p xahcnormal
#SBATCH -o log/%j
#SBATCH -e log/%j

export core_n=32
export filename=freecad
export fileformat=stp
export numlevels=1
export numrefine=1
export maxh=1000.0 
export minh=0.0
export input_path=./$filename\.$fileformat
export output_path=../output8/$filename\_r$numrefine\_l$numlevels\_$core_n\_max$maxh\_min$minh\/

mpicc --version
mpirun --version

module list 

mpirun -disable-auto-cleanup -np $core_n ./mesh_occ_mpi -i $input_path -o $output_path -l $numlevels -r $numrefine --maxh $maxh --minh $minh -v -adj