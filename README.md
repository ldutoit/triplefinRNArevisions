# README

This project is associated to:

Riordan, B., King, T., Dutoit, L., Gemmell, N., Hickey, A., & Johnson, S. Genetic and mitochondrial contributions to thermal resilience in an estuarine and intertidal fish. Unpublished.

Feel free to ask for help or any additional file to breana.marie2@gmail.com when accessing this repository.

# Summary

We built the transcriptome of the estuarine triplefin (_Forsterygion_nigripenne_) and the common triplefin (_Forsterygion_lapillum_) denovo using Trinity. We then performed gene annotation using Trinotate and further downstream analyses in R Studio v4.3.0). Using this repository, the analysis should be entirely reproducible. The main output files are here.

# revisions

Two new transcriptomes, one for nigripenne and one for lapillum were assembled with [Trinity_denovo_transcriptome_species_specific.md](Trinity_denovo_transcriptome_species_specific.md).

Raw reads for obtaining matrix counts were reprocessed through [Quality_control_and_trimming.md](Quality_control_and_trimming.md) followed by [Obtaining_iscount_matrix.md](Obtaining_iscount_matrix.md). That last file was run three times, one on the hybrid assembly from the original submission, one on the species-specific nigripenne assembly and one on the species-specific lapillum assembly.

The transcriptome assemblies are two big to be uploaded here, but they are at /nesi/project/uoo00116/Breana/.

The gene counts are:

[RSEM_gene_counts_hybrid.txt](RSEM_gene_counts_hybrid.txt): gene counts from scratch on Breana's transcriptome  assembly
[RSEM_gene_counts.txt](RSEM_gene_counts.txt): Breanna's original submission gene counts from MolEcol submission (found on github)
[RSEM_gene_counts_L.txt](RSEM_gene_counts_L.txt): gene counts from scratch on lapillum transcriptome  assembly
[RSEM_gene_counts_N.txt](RSEM_gene_counts_N.txt): gene counts from scratch on nigripenne transcriptome  assembly

I also ran OrthoFinder in [orthofinder.md](orthofinder.md) between two or three assemblies to be able to do species specific assemblies.

Finally,DE Analysis for the brain of lapillum only, but on all three transcriptomes (original hybrid transcriptome  one with both gene counts matrix) is in [DE_Analysis_brain_lapillum](DE_Analysis_brain_lapillum).
