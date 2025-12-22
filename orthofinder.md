## transdecoder to get proteins...

```
/home/ludovic.dutoit/revisions/comparing_assemblies
```

```
#!/bin/sh
#SBATCH --job-name=transdecoder
#SBATCH --account=uoo03946
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --hint=nomultithread
module load TransDecoder
TransDecoder.LongOrfs -t RF_trinity_outputL.Trinity.fasta                    
TransDecoder.LongOrfs -t RF_trinity_outputN.Trinity.fasta                    
TransDecoder.LongOrfs -t RF_trinity_output.Trinity.fasta  
TransDecoder.Predict -t RF_trinity_outputL.Trinity.fasta
TransDecoder.Predict -t RF_trinity_outputN.Trinity.fasta
TransDecoder.Predict -t RF_trinity_output.Trinity.fasta
```
## into orthofinder

```
module load OrthoFinder
mkdir orthofinder_run
cp *pep orthofinder_run

```

```
#!/bin/sh
#SBATCH --job-name=orthofinder
#SBATCH --account=uoo03946
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --hint=nomultithread
module load OrthoFinder
orthofinder -f orthofinder_run/
```

the results we want is what is in OrthoGroups/Orthogroups_single_orthologous.., those have all the groups where we have one group per species.
