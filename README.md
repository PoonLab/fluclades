# Genetic taxonomy of influenza A virus subtypes

Scripts and output files associated with a manuscript "Prospects for a sequence-based taxonomy of influenza A virus subtypes".

Input data files can be obtained from Zenodo under a Creative Commons license:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8119571.svg)](https://doi.org/10.5281/zenodo.8119571)

## Scripts

* `chainsaw-plot.R` - R script to visualize the number of subtrees produced by the edgewise clustering ("chainsaw") method as a function of the internal branch length cutoff.  
This script is applied to results obtained for protein sequence phylogenies reconstructed for all eight influenza A virus (IAV) genome segments, with an emphasis on hemagglutinin (HA) and neuraminidase (NA) proteins.


* `chainsaw.py` - Python script implementing the edgewise clustering method.  Requires Biopython.
Running the script without any arguments prints a histogram summary of branch lengths to the console.
Specifying a branch length cutoff with `--cutoff` prints a summary of the resulting subtres (defaulting to `-f summary`).
Setting the `-f` option to `labels` writes a detailed CSV output listing subtree assignments for all tips.
The script also calculates the normalized mutual information between the subtree partition and subtype labels.

* `compress-seqs.py` - This Python script looks for exact matches in unaligned sequences of the input FASTA file, and writes the unique sequence to an output FASTA file using the first label encountered.  All other duplicate labels are written to a CSV file to link them to the first label.  This script also filters out sequences with an excessive number of ambiguous amino acids (`X`).

* `concat-genes.py` - This Python script concatenates the non-overlapping amino acid sequences for M1/M2 or NS1/NS2 records from the same isolate.  The input is assumed to be a CDS FASTA file generated by the NCBI Genbank interface.

* `filter-prot.py` - This Python script applies an initial filter on the CDS FASTA files downloaded from Genbank.  It uses regular expressions to remove records that do not correspond to the query protein.

* `get-metadata.py` - The default sequence names for Genbank CDS downloads are not very informative, so this script is used to retrieve more useful metadata such as the strain name and collection date from the database based on the accession number.  It takes either a FASTA or NWK file as input.  The results are written to a CSV file.

* `midpoint.R` - This small R script simpily calls the `midpoint` rooting function of the [`phangorn`](https://cran.r-project.org/web/packages/phangorn/index.html) package on the input tree.

* `plot-trees.R` - This R script requires the R package [`ggfree`](https://github.com/ArtPoon/ggfree).  It generates plots of the large HA and NA phylogenies, colouring branches based on subtype labels on the tips.

* `relabel-fasta.py` - This Python script uses the CSV generated by `get-metadata.py` to replace the sequence names in the user-specified FASTA input file.

* `subtree-grid.R` - This R script was used to generate the supplementary figure summarizing the results of node-wise clustering of the HA phylogeny.

* `subtyping.py` - This Python script implements the nodewise clustering method, calculating a number of summary statistics for every internal node of the input tree.


## Results

