#!/usr/bin/zsh

### Job name
#SBATCH --job-name=PCMLatticeSolver

# declare the merged STDOUT/STDERR file
#SBATCH --output=outputsolver.%J.txt

### Request the time you need for execution. The full format is D-HH:MM:SS
### You must at least specify minutes or days and hours and may add or
### leave out any other parameters
#SBATCH --time=00:20:00
  
#### Request the memory you need for your job. You can specify this
### in either MB (1024M) or GB (4G).
#SBATCH --mem-per-cpu=250M

### Load the required module files
module load GCC/10.3.0
module load OpenMPI/4.1.1 
module load OpenFOAM/9

# Execute the wmake command
wmake
