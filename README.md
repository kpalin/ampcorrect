# ampCorrect software for correcting and collapsing nanopore sequenced amplicons

This repository contains two scripts assisting working on nanopore
sequenced amplicons. First is `polish_amplicon.sh` which takes fasta
file of nanopore reads, clusters them according to similarity,
corrects by multiple alignment and maps the alignment consesuses to
reference genome.

This protocol is heavily influenced by [nanocorrect](https://github.com/jts/nanocorrect)

	usage:
	./polish_amplicon.sh [-c WORKDIR] [-r REF.fasta] input.fasta

	-c WORKDIR   Continue work on this directory
	-r REF.fasta Genome reference for mapping the final results. 
	-h           Show this message and exit.



Another script is used for separating amplicon reads that are almost
equal and are separated by discovery of short "splitter" reads.

	usage:
	./split_reads.sh -F splitters.fasta sequence.(fast5|fasta)

	Map short sequences in splitters.fasta to reads in sequence.fasta and
	generate SAM like output file usable for e.g. split_fasta_by_sam.py
	
	-F splitters.fast(a|5)  Split the reads according to alignment to these.
	-d          Use default lastal arguments, instead of ones for nanopore reads
	-h          Show this message and exit.


These scripts have dependency on
[sumaclust](https://git.metabarcoding.org/obitools/sumaclust/wikis/home)
and [poaV2](https://sourceforge.net/projects/poamsa/). These need to
be in $PATH.

If you use this software in academic publications, please cite

> “Detection and Analysis of Somatic LINE-1 Insertions in Colorectal Cancer by Long-distance Inverse-PCR and Nanopore Sequencing” by Pradhan et al. (Submitted)
