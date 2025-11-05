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
run trinitynfrun.sh

```
#!/bin/sh
module load Nextflow fastp FastQC Trinity Java/17
nextflow run nf-core/denovotranscript --input samples_filenfL.csv --outdir RF_trininity_nfL -r main --profile myslurm
```

```
 sbatch -A uoo00116 -t 2-00:00:00 -c 1 --mem 1G trinitynfrun.sh
```
