#!/bin/bash                                                                           
## The above line indicates that script should be interpreted using                   
## Bash shell                                                                         

## The following script is used to update the simulation from the previous saved file\
                                                                              

## Slurm Directives ##                                                                

#SBATCH -J Looping_angle                ##Short for --jobname                         
#SBATCH -N 1                ##Short for --nodes                                       
#SBATCH -n 56               ##Short for --nstasks #no of tasks                        
#SBATCH -t 02:30:00         ##Short for --time #wall time limit                       
#SBATCH -o slurm-%j.out     ##Short for --output #Standard output name                
#SBATCH -e slurm-%j.err     ##Short for --error  #Standard error name                 
#SBATCH -p condo            ##Short for --partition #partition name                   
#SBATCH -q condo            ##Short for QOS name                                      
#SBATCH -A csd837           ## Allocation name                                        
#SBATCH --mail-type END     #Optional, Send mail when job ends                        
#SBATCH --mail-user <email> #Optional, send email to this address                     
#SBATCH --mem=128G

## Execution Block ##                                                                 

mpirun -np 56 ./main_exe

