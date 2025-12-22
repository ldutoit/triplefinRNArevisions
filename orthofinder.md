## predicting orf

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

then orthofinder somehow
