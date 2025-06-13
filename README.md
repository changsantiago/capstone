# Computational and Systems Biology Capstone Research
Undergraduate capstone research in computational and systems biology, hosted by the Junction of Statistics and Biology Lab at UCLA.

## Source files:
Genome sequences (FASTA) and annotation features (GTF):\
NCBI RefSeq assembly and annotation for WBcel235 found at this [site](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002985.6/). Download [here](https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000002985.6/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GTF).

phyloP scores (BigWig):\
UCSC Genome Browser basewise conservation scores (phyloP) of 134 nematode genomes with C. elegans for WBcel235 found at this [site](https://hgdownload.soe.ucsc.edu/goldenPath/ce11/phyloP135way/). Download [here](https://hgdownload.soe.ucsc.edu/goldenPath/ce11/phyloP135way/ce11.phyloP135way.bw).

## Code:
script.py: extracting protein-coding mRNA strands with a 5' UTR from source files\
script2.py: constructing data frame with position and strand for each phyloP score\
script3.py: identifying and removing polycistronic genes\
script4.py: identifying and recording count and indices of type A, B, C uORFs\
script5.py: extracting cDNA sequence, distance from CDS, and Kozak sequence for uORFs\
script6.py: calculating nucleotide frequency by position based on CDS and calculating Kozak similarity scores for uORFs\

## Data:
nonpolycistronic.tsv: non-polycistronic protein-coding genes with a 5' UTR\
uORFs.tsv: type A uORFs\
closefaruORFs.tsv: close and far type A uORFs
