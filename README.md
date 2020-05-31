# Linear Mixed Model for single trait eQTL discovery

<br />

This is a single trait and easy to run implementation of the linear mixed model [LIMIX](https://www.biorxiv.org/content/10.1101/003905v2) framework for QTL mapping. LIMIX is an efficient python implementation of [Linear Mixed Models](https://en.wikipedia.org/wiki/Mixed_model) used to properly control for [population structure](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007309) in association studies.

This single trait implementation can perform univariate association tests between a single phenotype and multiple variants. A major feature is the ability to compute empirical p-values based on phenotype permutations.
Empirical p-values are the most reliable way to assess association because they make no assumption on the genotype and phenotype distributions.

To improve parallelization, this script performs test on a single gene (that has to be provided in the options).


<br /><br />


## Prerequisites

<br />

Install [Anaconda 2](https://www.anaconda.com/products/individual) following the instruction on the Anaconda website. LIMIX is based on python2.


<br /><br />

## Installation

<br />

Clone this GitHub repository 

<pre><code>    git clone https://github.com/mforneris/LIMIX-single-trait.git
</code></pre>


Enter in the LIMIX implementation folder

<pre><code>   cd LIMIX-single-trait
</code></pre>

Create an Anaconda environment that includes all the packages necessary to run LIMIX

<pre><code>    conda env create -f config/environment.yml -n LIMIX_env
</code></pre>

Activate the conda environment

<pre><code>    source activate LIMIX_env
</code></pre>

Finally install vcftools

<pre><code>    source install vcftools
</code></pre>


And you are good to go!


<br /><br />


## Expected inputs



<br />

#### 1. Genotypes

The genotype file contains the variants information that will be used as independent variable in the test. Most variant files are in [vcf format](https://samtools.github.io/hts-specs/VCFv4.2.pdf). This LIMIX uses [hdf5 format](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) files as inputs to make processing faster and data storage more compact. It is adviced to filter variants with maf < 5% and missing values > 20%.

For a guide on how to convert vcf files to hdf5 look at the command line section.

<br />

#### 2. Phenotypes

The phenotype file contains the quantitative trait information that will be used as dependent variable in the test. The phenotype matrix should have individuals as rows and phenotypes as columns. The single triat LIMIX can perform tests for each phenotype singularly. Gene expression values can be corrected for confounding factors using [PEER](https://www.nature.com/articles/nprot.2011.457).

The command line section contains a guide on how to convert csv files to hdf5.

<br />

#### 3. Neutrally evolving variants to infer population structure

Neutrally evolving variants are used to infer population structure. It is important to distinguish between variants under selective pressure and variants that can be used as a proxy of time. Neutrally evolving variants can be found as coding synonnymous variants or variants in short introns (>60bp).

The vcf including neutrally evolving variants shloud be converted in hdf5 format.

<br />

#### 4. Range of variants indexes for each gene

Thsi file should contain three TAB-separated columns. The first column is the gene name, the second and the third should contain the indexes of the first and last variant to be tested for association with the gene. The variant indexes start from 0 (with the first variant in the vcf). The interval should be left-closed, right-open [a,b).

E.g. if I want to test variants from the 500 to the 673 for association with gene ENSG0000, the corresponding entry in the file will be:

ENSG0000  &nbsp; &nbsp; &nbsp; &nbsp; 499  &nbsp; &nbsp; &nbsp; &nbsp; 673

<br /><br />


## Output 


The oput file will contain 6 columns:
1. The chromosome of the variants
2. The position (1-based) of the variant
3. The gene name
4. The p-value of the linear mixed model test
5. The beta of the linear mixed model test
6. The empirical p-value

<br />



## Command line examples


<br />

#### Convert vcf files to hdf

First you need to convert the vcf file in 012 format using vcftools

<pre><code> vcftools \
    --gzvcf example/input_data/test_variants.vcf.gz \
    --012 --out example/input_data/test_variants
</code></pre>

Then you can use limix_converter to convert the 012 file to hdf

<pre><code> python2.7 bin/limix_converter.py \
    --g012=example/input_data/test_variants.012 \
    --outfile=example/input_data/test_variants.hdf
</code></pre>

The same procedure can be applied to convert the population structure vcf file.

<br />

#### Convert csv phenotype file files to hdf

To convert the gene expression file into hdf format use:

<pre><code> python2.7 bin/limix_converter.py \
    --csv=example/input_data/gene_expression.csv  \
    --outfile=example/input_data/gene_expression.hdf
</code></pre>

Make sure that the expression file is properly formatted.


<br />

#### Run LIMIX test


Here is an example on how to run limix single trait test on a gene. The gene name is specified with the option -x. The number of permutations to compute the empirical p-value is specified with the option -n. The number of permutations will influence the quality of the empirical p-value estimate. 

<pre><code> python2.7 bin/limix_single_trait_gene_based.py \
        -g example/input_data/test_variants.hdf \
        -p example/input_data/gene_expression.hfd5 \
        -r example/input_data/population_structure_variants.hdf \
        -c example/input_data/variants_indexes_to_genes.txt \
        -x FBgn0000024 \
        -n 10000 \
        > example/output_data/limix_single_trait_eQTLs_FBgn0000024.txt
</code></pre>

You can find a command line example in the bin folder.




