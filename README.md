# cghla

A set of tools that are useful for finding the MHC/HLA alleles for a given
sample using WGS, WES, or RNAseq data. In particular, these tools are designed for
the identification of MHC alleles to use in identification of novel peptides that can
be presented by MHC Class I proteins for finding presented neoepitopes in cancer.
As such, the specificity of the results are designed only for this purpose.

These tools implement a particular workflow that assumes you have the following:

    bwa
    samtools
    ngsutilsj

Additionally, the tools are built to work with data from the IPD-IMGT/HLA database
which is available from: https://hla.alleles.org

Finally, final prediction can use MHC peptide binding affinities from mhcFlurry which
are available here: https://openvax.github.io/mhcflurry-motifs/mhcflurry.allele_sequences.csv

The basic workflow can be described as:

## WGS
1) Find the location of HLA alleles in the reference genome.
2) Extract flanking sequence from the reference genome, and adding this sequence
   to the allele sequences
3) Extraction of mapped reads to the HLA locus in an already aligned BAM file
4) Re-alignment of these reads (optionally with unmapped reads as well) to the 
   now extended flanking HLA allele sequences
5) Score HLA allele pairs (maternal/paternal) to find the most likely genotypes
   for the individual.
6) Predict the most likely genotypes

## RNAseq
The RNAseq/WES workflow is similar, but instead of using the `*_gen.fasta` versions
of the IPD-IMGTR/HLA sequences you can use the `*_nuc.fasta` versions. When using these
coding mRNA sequence, you can skip steps 1 and 2. And instead of extracting a sub-set
of reads, you can align the raw FASTQ sequence. The scoring and prediction steps remain
the same.

