## Introduction
The following repository contains the scripts used to process that data as well
as generate the figure in "Leveraging histone modifications to improve genome annotations" - Mendieta et al 2021. 

The corresponding directories are broken down by script type, with the `snakemake` directory being further broken down by analysis type. 

## Annotating Genomes Using Chromatin Modification Data (Script)
Below is a detailed example of how to utilize the annotation script implemented
in Mendieta et al. Detailed infomration about pre-processing ChIP-seq data to
utilize this script is included in the script `update_annotation_by_chromatin_mods.snake`. 
But paired down, the script requires a set a genome annotations in bed format
(genes and lncRNAs). A set of peaks from ChIP-seq peak calling focused on modifications 
corresponding to elongation which should occur throughout the gene body
(H3K36me3, H3K4me1), as well as initiation modifications which correlate with
the promoter of a gene (H3K4me3, H3K56ac). At least one modification of
each type is required. Additionally stranded RNA-seq data is required to ensure
that non-expressed gene features overlapping our modifications will be ignored.
Finally a bigwig file corresponding to one of the elongation mods is helpful,
as it allows the script to identify large genes while accounting for potential
mappability issues, as well as lack of upstream signal from the elongation mod.

Basic Script Command: 
```
python scripts/Update_annotation.py -broad {input.broad} \
-narrow {input.narrow_peaks} -annotation {input.gene_annotation} \
-bw {input.H3K36me3_bw} \
-RNA {input.rna_file} -lncRNA {input.gene_annotation} -o {output}
```

A list of test/example files can be found in this github repo under the file
`example_data`.

The output of this script is a list of annotations - with their potential
annotation corrected class append as original.annot_ANNOTATION_CLASS. Potential
novel annotations are awarded the name `Possible_novel_gene`. 

Note: This script is implemented in python 3.6, and requires pybedtools, and
pysam. Both of which can be easily installed using either `homebrew` or
`Anaconda3`. Additionally since many genomes which this analysis can be run on
are fairly large, they will need to be broken up on a chromosome by chromosome
basis. This ensures that pysam is able to read the full RNA-seq BAM file
without overflowing. 


## Annotating Genomes Using Chromating Modification Data (Pipeline)

The implemenation of this method is further detailed in the snakemake script
located `snakemake_scripts/12.Annotate_other_species`. Here you will find an
example implementation in soybean (glycine max). Structure of the inputs can be 
gleaned from the snakemake config `yaml` file,  as well as a `sample` 
file corresponding to the features mentioned above.
Note that tissues need to be matched in the `sample` file. This `snakemake`
analysis was done in version of `5.7.1`. Since then aggregate command behavior
is slighly modified, so be aware that running this analysis on more updated
version of snakemake may result in unexpected behaviors.

