The following code is to be utilised for shorter reads (2x 51 bp paired-end reads). For longer pair-end reads to be used for generating a denovo transcriptome, see the file [QC_and_trimming_long](https://github.com/breanariordan/triplefinRNA/blob/main/QC_and_trimming_long.md)

# Moving raw files

```
pwd                       # past working directory, find where you are
cd OG7591-691868198       # identify/move to to directory where raw RNASeq files are
cd ../nobackup            # move back to nobackup folder
mkdir sourcefiles         # make a new directory called 'sourcefiles'
cd sourcefiles
mkdir run                 # make a new directory for the raw reads
mkdir re-run
cd ..                     # move back to nobackup directory
tar xvf OG7591-39750650.tar

cp Desktop/OG7591-397506150/OG7591-691868198/*/*gz sourcefiles/run/
cp Desktop/OG7591-397506150/OG7591_re-run-694316623/*/*gz sourcefiles/rerun/
```

# Quality control

```
cd sourcefiles                # if you aren't already in it
mkdir QC                      # make a folder called QC
cd run
module load FastQC/0.12.1     # load fastqc for quality control
fastqc -o QC/ run/*           # perform quality control on the raw data
cd ../QC                      # move into the QC directory
for filename in *.zip         # unzip fastqc text file outputs
    do
    unzip $filename
    done
cat */summary.txt > ~/nobackup.sourcefiles/QC/fastqc_summaries.txt # generate summary of QC stat tests to see where failures have occurred
module load MultiQC/1.13-gimkl-2022a-Python-3.10.5 # load MultiQC
cd ..                         # move back to sourcefiles
mkdir MultiQC                 # make a directory called 'MultiQC'
cd MultiQC
cp ../QC/* ./                 # move all QC files for individuals samples into MultiQC directory
multiqc .                     # perform multiqc analysis 
```

# Trimming
if fastqc and multiqc results indicate it is needed
this will use the cutadapt programme to trim low quality reads out of the files and put them into a new folder called 'Trimmed'

```
cd /home/ludovic.dutoit/revisions/OG7591-691868198/sourcefiles
mkdir Trimmed                 # make directory 'Trimmed'
module load cutadapt/4.1-gimkl-2022a-Python-3.10.5    # load latest cutadapt version
for filename in run/*fastq.gz
    do
    base=$(basename ${filename} .fastq.gz)
    cutadapt -q 20 -a AGATCGGAAGAG -o Trimmed/${base}.trimmed.fastq.gz ${filename} > Trimmed/${base}.log
    done
for filename in rerun/*fastq.gz
    do
    base=$(basename ${filename} .fastq.gz)
    cutadapt -j 2 -q 20 -a AGATCGGAAGAG -o Trimmed/${base}.trimmed.fastq.gz ${filename} > Trimmed/${base}.log
    done

cd ../MultiQC
cp ../Trimmed/*log .          # take all the trimmed .log files
multiqc .                     # and run the multiqc again
```
```
# /home/ludovic.dutoit/revisions/OG7591-691868198/sourcefiles/Trimmed

#concatenate them
	mkdir -p concatenated_fastq

	for s in $(seq 1 52); do
	    s_tag="S${s}"
	    
	    files=$(ls *_${s_tag}_L001_R1_001.trimmed.fastq.gz 2>/dev/null)
	    
	    if [ $(echo "$files" | wc -l) -eq 2 ]; then
	        out="concatenated_fastq/${s_tag}_L001_R1.trimmed.fastq.gz"
	        zcat $files > "$out"
	        echo "Concatenated $s_tag"
	    else
	        echo "Skipping $s_tag (found $(echo "$files" | wc -l) files)"
	    fi
	done
```

Renaming to be done

Stuff to compare species specific mapping rate to each transcriptome.
