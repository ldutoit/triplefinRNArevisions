## compring_run_rerunonly.md

We want to check if the original analyses was run on the big rerun only or the concatenated files

First compare concatenated:
```
less ~/revisions/comparing_assemblies/10N2heart/RSEM.genes.results | head -n 2
#TRINITY_DN0_c0_g1       TRINITY_DN0_c0_g1_i1,TRINITY_DN0_c0_g1_i2,TRINITY_DN0_c0_g1_i3,TRINITY_DN0_c0_g1_i4,TRINITY_DN0_c0_g1_i5,TRINITY_DN0_c0_g1_i6   216.97  167.14  43.84 8
.61     8.68
```

to ms results:
```
less /nesi/nobackup/uoo03946/denovo_transcriptome/RSEM/10N2heart.txt	
#TRINITY_DN0_c0_g1       TRINITY_DN0_c0_g1_i1,TRINITY_DN0_c0_g1_i2,TRINITY_DN0_c0_g1_i3,TRINITY_DN0_c0_g1_i4,TRINITY_DN0_c0_g1_i5,TRINITY_DN0_c0_g1_i6   217.03  167.11  44.87 8
.83     8.88
```
Almost the same.

What about rerun only :

We'll go take 10 N2 heart but just the rerun it is S31

```
mkdir checkfullorrerun
cd  checkfullorrerun
cp /home/ludovic.dutoit/revisions/OG7591-691868198/sourcefiles/Trimmed/AAC7KWHHV-7591-31-49-1_S31_L001_R1_001.trimmed.fastq.gz 10_N2_heart.fastq
cp ../comparing_assemblies/RF_trinity_output.Trinity.fasta* .
```

make a samples_file_small.txt:

```10N2heart  10N2heart  10_N2_heart.fastq```

```
cp ../comparing_assemblies/iscountmatrix_hybrid.sl .
sbatch iscountmatrix_hybrid.sl
```

it is 3788765
