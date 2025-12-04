
### Comparing assemblies

```
#cd /home/ludovic.dutoit/revisions
mkdir comparing_assemblies
cd comparing_assemblies
cp /nesi/nobackup/uoo03946/revisions/RF_trinity_outputN.Trinity.fasta . # N
cp /nesi/nobackup/uoo03946/revisions/RF_trinity_outputL.Trinity.fasta . # L
cp /nesi/project/uoo03946/denovo_transcriptome/RF_trinity_output.Trinity.fasta . #hybrid
```

## rnaQUAST


```
module load rnaQUAST
rnaQUAST.py --transcripts RF_trinity_output.Trinity.fasta --o hybrid # jus testes not ran
```
...
