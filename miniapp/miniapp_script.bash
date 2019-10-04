### LSF syntax
#BSUB -nnodes 1                   #number of nodes
#BSUB -W 60                      #walltime in minutes
#BSUB -G guests                   #account
#BSUB -e exa_error.txt              #stderr
#BSUB -o exa_out.txt            #stdout
#BSUB -J ExaCMech_Miniapp               #name of job
#BSUB -q pdebug                   #queue to use

#The OpenMP options could probably be played around with
#The first option just checks the weak scaling of the code
#if we do 7x the work and have 7x the cores available do
#we run in the same amount of time compared to the serial job
#On Summit, we have 42 cores available each with 4 threads and
#6 GPUs available. The current recommendations for splitting up work
#is to have one MPI rank per GPU (1 MPI rank corresponds to 7 cores
#if we make use of OpenMP)
export OMP_NUM_THREADS=7
export OMP_PROC_BIND=close
export OMP_PLACES=cores
#We want to make sure only 1 gpu is actually being used.
export CUDA_VISIBLE_DEVICES=0
./orientation_evolution ./option_cpu.txt
./orientation_evolution ./option_openmp.txt
#Here we're checking the strong scaling where if we increase
#the number of threads by 4 times is our runtime 1/4 of the original.
#If it isn't there's a good chance we're limited by the cache memory
#that each thread shares for an individual core.
export OMP_NUM_THREADS=28
export OMP_PROC_BIND=close
export OMP_PLACES=threads
./orientation_evolution ./option_openmp.txt
#Now we can see just what sort of scaling we get on the GPU
./orientation_evolution ./option_cuda.txt
