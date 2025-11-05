## NextFlow pipeline

It runs here:

```
cd /nesi/nobackup/uoo03946/revisions
```
with the following nextflow.config:

```
profiles {
    myslurm {
        profile {
            // Use SLURM executor
            process.executor = 'slurm'
            // Default CPUs and memory per process
            process.account = 'uoo00116'
        }
    }
}
```
and the batch file trinitynfrun.sh.

```
#!/bin/sh
module purge
module load SAMtools Nextflow fastp FastQC Trinity Java/17 
nextflow run nf-core/denovotranscript --input samples_filenfL.csv --outdir RF_trinity_nfL  --assemblers trinity -r main -profile myslurm -resume
```

Run the pipeline using:
```
sbatch -A uoo00116 -t 2-00:00:00 -c 24 --mem 72G trinitynfrun.sh
```
