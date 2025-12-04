
### Comparing assemblies

```
#cd /home/ludovic.dutoit/revisions
mkdir comparing_assemblies
cd comparing_assemblies
cp /nesi/nobackup/uoo03946/revisions/RF_trinity_outputN.Trinity.fasta . # N
cp /nesi/nobackup/uoo03946/revisions/RF_trinity_outputL.Trinity.fasta . # L
cp /nesi/project/uoo03946/denovo_transcriptome/RF_trinity_output.Trinity.fasta . #hybrid
```

remap both species to each assemblies...

/Obtaining_iscount_matrix.md


## rnaQUAST


```
module load rnaQUAST
rnaQUAST.py --transcripts RF_trinity_output.Trinity.fasta --o hybrid # jus testes not ran
```
...


I’ve only used this once, but it worked well.  It can do the same first option in transrate,  remaping raw reads and detect artefacts. If you have a reference genome, it can also report missing, mapped, and duplicated genes base on the genome. Although you can run the two transcriptomes at the same time. It does not have the option of map transcriptome vs transcriptome like Transrate. But, this tool is more recent and better maintained.

 
