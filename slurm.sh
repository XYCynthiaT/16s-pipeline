#!/bin/bash
#SBATCH --job-name=16s-pipeline 
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --time=1-0 
#SBATCH --mem=2000 
#SBATCH --partition=high 
#SBATCH --output=slurmout/snakemake_%A.out 
#SBATCH --error=slurmout/snakemake_%A.err 
##SBATCH --mail-type=ALL
##SBATCH --mail-user=myemail@email.com

snakemake -j 20 --cluster-config config/cluster.json --cluster "sbatch -J {cluster.name} -p {cluster.partition} -N {cluster.nodes} -n {cluster.ntasks} --mem={cluster.mem} -t {cluster.time} -o {cluster.output} -e {cluster.error}"