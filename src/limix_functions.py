from limix_load_packages import *

def read_phenotypes(hdf_file):
        hdf_phenotype = h5py.File(hdf_file, 'r', driver='core')
        phenotypes_matrix = hdf_phenotype[u'phenotype'][u'matrix'][()]
        samples_ID=list(hdf_phenotype[u'phenotype'][u'row_header'][u'sample_ID'])
        phenotypes_ID = list(hdf_phenotype[u'phenotype'][u'col_header'][u'phenotype_ID'])
        phenotypes = pd.DataFrame(data = phenotypes_matrix, index = samples_ID, columns = phenotypes_ID)
        return phenotypes

def impute_genotypes(genotypes_matrix):
        """
        This function imputes the missing genotypes. It assumes that missing genotypes have a value of 255 (or -1).
        The first step computes the average value of the genotyped samples. The second substitutes the imputed to the missing values.
        """
        genotypes_matrix_float = genotypes_matrix.astype(float)
        for index in range(genotypes_matrix_float.shape[1]):
                genotyped_lines_average = np.extract(genotypes_matrix_float[:,index]!=255.0, genotypes_matrix_float[:,index]).mean()
                np.place(genotypes_matrix_float[:,index], genotypes_matrix_float[:,index]==255.0, genotyped_lines_average)
        return genotypes_matrix_float

def read_genotypes(hdf_file):
        hdf_genotype = h5py.File(hdf_file, 'r', driver='core')
#        chrom = hdf_genotype[u'genotype'][u'col_header'][u'chrom'][()].tolist()
#        pos = hdf_genotype[u'genotype'][u'col_header'][u'pos'][()].tolist()
	pos_file_name = re.sub('.hdf[5]?', r'.012.pos', hdf_file)
	try:
		text_file = open(pos_file_name, "r")
	except:
		print "\nERROR!\nFile '%s' does not exist.\nThe script had to circumvent a bug from LIMIX and needs to read the position name directly from the 012 file.\nCheck the function file\n\n" %(pos_file_name)
	lines = text_file.readlines()
	chrom = [l.split('\t')[0] for l in lines]
	pos = [l.split('\t')[1].split('\n')[0] for l in lines]
	genotypes_matrix = hdf_genotype[u'genotype'][u'matrix'][()]	
        genotypes_matrix_imputed = impute_genotypes(genotypes_matrix)
        samples_ID = hdf_genotype[u'genotype'][u'row_header'][u'sample_ID'][()]
        genotypes_ID = [chrom[i] + "_" + str(pos[i])  for i in range(0, len(list(chrom)))]
        genotypes = pd.DataFrame(data = genotypes_matrix_imputed, index = samples_ID, columns = genotypes_ID)
        return genotypes

def read_cis(cisfile):
        cis_dict = {}
        with open(cisfile, 'r') as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                        val = [row[1], row[2]]
                        key = row[0]
                        cis_dict[key] = val
        return cis_dict

def compute_relatedness_matrix(genotypes):
        Xc = np.array(genotypes)
        #1. standardization
        Xc -= Xc.mean(0)
        Xc /= Xc.std(0)
        #2. dot product
        Kinship = sp.dot(Xc, Xc.T)
        #3. normalization
        Kinship /= Kinship.diagonal().mean()
        return Kinship

def extract_variants(genotypes, cis_snp_dict, gene):
	pattern = re.compile("FBgn")
	g = [i for i in gene.split("_") if pattern.match(i)][0]
	idx_start = int(cis_snp_dict[g][0])
	idx_end = int(cis_snp_dict[g][1])
	genotype_subset = genotypes.iloc[:, idx_start:idx_end]
	return genotype_subset

def monovariate_LMM(genotypes_subset, gene_expression, relatedness_matrix):
	monovariate_lmm = qtl.qtl_test_lmm(genotypes_subset.values, gene_expression.values, K=relatedness_matrix)
	betas = monovariate_lmm.getBetaSNP()
	pvalues = monovariate_lmm.pvalues
	return(monovariate_lmm, betas, pvalues)

def compute_permuted_pvalues(genotypes_subset, gene_expression, relatedness_matrix, n_permutations):
	permuted_pvalues_monovariate = []
	for i in range(n_permutations):
		sys.stderr.write("      ... permutation %d\n" % i)
		phenotype_index = range(len(gene_expression))
		shuffle(phenotype_index)
		shuffled_gene_expression = gene_expression.iloc[np.asarray(phenotype_index)]
		shuffled_gene_expression.index = gene_expression.index
		permuted_monovariate_lmm, permuted_betas, permuted_pvalues = monovariate_LMM(genotypes_subset, shuffled_gene_expression, relatedness_matrix)
		min_pvalue_in_permutation = permuted_pvalues.min()
		permuted_pvalues_monovariate.append(min_pvalue_in_permutation)
	return permuted_pvalues_monovariate

def empirical_pvalues(n_permutations, permuted_pvalues_monovariate, pvalues):
	empirical_single_effect_pvalues = sp.zeros((1, pvalues.shape[1]))
	permuted_pvalues = np.array(permuted_pvalues_monovariate)
	for s in range(pvalues.shape[1]):
		pvalue = pvalues[:,s]
		permuted_pv_lower_than_pvalue = np.extract(permuted_pvalues < pvalue, permuted_pvalues).shape[0]
		empirical_pvalue = float((permuted_pv_lower_than_pvalue + 1)) / (n_permutations + 1)
		empirical_single_effect_pvalues[0, s] = empirical_pvalue
	return empirical_single_effect_pvalues




