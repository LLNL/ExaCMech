### LSF syntax
#BSUB -nnodes 1                   #number of nodes
#BSUB -W 60                      #walltime in minutes
#BSUB -G guests                   #account
#BSUB -e exa_error.txt              #stderr
#BSUB -o exa_out.txt            #stdout
#BSUB -J ExaCMech_Miniapp               #name of job
#BSUB -q pdebug                   #queue to use

#The OpenMP options could probably be played around with
export OMP_NUM_THREADS=28
export OMP_PROC_BIND=close
export OMP_PLACES=threads
#We want to make sure only 1 gpu is actually being used.
export CUDA_VISIBLE_DEVICES=0
./orientation_evolution ./option_cpu.txt
./orientation_evolution ./option_openmp.txt
./orientation_evolution ./option_cuda.txt
