#!/usr/bin/env python

from limix_load_packages import *
from limix_functions import *

# Setting up arguments
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--geno_file', dest='genoin', help="Input genotype files must be in hdf5 format")
parser.add_argument('-p', '--pheno_file', dest='phenoin', help="Input phenotype file must be in hdf5 format")
parser.add_argument('-r', '--relatedness_hdf', dest='relatedness_hdf', default=None, help="Input genotype files contaning only neutral variants. This hdf file is used to infer the relatedness matrix. Be carefull not to use the entire set of variants since the population structure is hidden by selective pressure on variants.")
parser.add_argument('-c', '--cis_snp_file', dest='cisfile', help="Input cis landscape for each gene. 3 Columns file with gene name, Index (0 based) for first variant associated to gene, index for last variant associated to gene")
parser.add_argument('-n', '--n_permutations', dest='permut', type=int, help="Number of permutations to perform to compute empyrical p-values")
parser.add_argument('-x', '--gene_to_test', dest='gene_to_test', default=None, help="The gene name on which to perform the test. This setting is used for better parallelization (one job for each gene and makes it easy to follow the outputs.")



def main():
	args = parser.parse_args()

	sys.stderr.write('... Reading inputs\n')
	phenotypes = read_phenotypes(args.phenoin)
	genotypes = read_genotypes(args.genoin)
	genotypes_relatedness = read_genotypes(args.relatedness_hdf)
	cis_snp_dict = read_cis(args.cisfile)
	gene = args.gene_to_test
	n_permutations = args.permut

	sys.stderr.write('... Input elaboration\n')
	relatedness_matrix = compute_relatedness_matrix(genotypes_relatedness)
	gene_names = list(phenotypes.axes[1])

	genotypes_subset = extract_variants(genotypes, cis_snp_dict, gene)

	sys.stderr.write('... Fitting model\n')
	gene_expression = phenotypes[gene]
	monovariate_lmm, betas, pvalues = monovariate_LMM(genotypes_subset, gene_expression, relatedness_matrix)

	sys.stderr.write("... Starting Permutations\n")
	permuted_pvalues_monovariate = compute_permuted_pvalues(genotypes_subset, gene_expression, relatedness_matrix, n_permutations)

	sys.stderr.write("... Computing Empyrical pvalues\n")
	empirical_monovariate_effect_pvalues = empirical_pvalues(n_permutations, permuted_pvalues_monovariate, pvalues)

	sys.stderr.write("... Printing output\n")
	print "%s\t%s\t%s\t%s\t%s\t%s" % ("Chromosome", "Position", "Gene", "Pvalue", "Beta", "Empirical_pvalue")
	chrom = [v.split('_', 1)[0] for v in list(genotypes_subset.columns.values)]
	pos = [v.split('_', 1)[1] for v in list(genotypes_subset.columns.values)]
	gene_name_list = [gene for v in range(genotypes_subset.shape[1])]
	pvalues_list = pvalues.tolist()[0]
	betas_list = betas.tolist()[0]
	empirical_pvalues_list = empirical_monovariate_effect_pvalues.tolist()[0]
	for i in range(len(chrom)):
		print "%s\t%s\t%s\t%.10e\t%.10e\t%.10e" % (chrom[i], pos[i], gene_name_list[i], pvalues_list[i], betas_list[i], empirical_pvalues_list[i])

main()



