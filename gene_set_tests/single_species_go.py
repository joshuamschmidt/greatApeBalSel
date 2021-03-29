import pysetperm as psp
import numpy as np
import pandas as pd
import pickle

n_perms = 100000
cores = 4
outpath='/Users/joshuaschmidt/Projects/greatApeBalSel/gene_set_tests/'
annotations = psp.AnnotationSet(annotation_file='../gene_annotations/hg19_nonRedundant_GOWINDA_noORS.bed', range_modification=2000)
tests = ['vips', 'kegg']
for test_set in tests:
    function_sets = psp.FunctionSets(function_set_file='../gene_sets/'+test_set+'_long_reformat.txt', min_set_size=5, annotation_obj=annotations)
    species = ['abelii', 'pygmaeus', 'gorilla', 'graueri', 'bonobo', 'troglodytes', 'schweinfurthii', 'ellioti', 'verus']
    for sp in species:
        candidates = psp.Variants(variant_file='../gene_set_test_files/'+sp+'_candidate_windows.txt')
        candidates.annotate_variants(annotation_obj=annotations)
        background = psp.Variants(variant_file='../gene_set_test_files/'+sp+'_background_windows.txt')
        background.annotate_variants(annotation_obj=annotations)
        test_obj = psp.TestObject(candidates,
                                  background,
                                  function_sets,
                                  n_cores=cores)
        permutations = psp.Permutation(test_obj, n_perms, cores)
        per_set = psp.SetPerPerm(permutations,
                                 function_sets,
                                 test_obj,
                                 cores)
        results = psp.make_results_table(test_obj, function_sets, per_set)
        results.to_csv(path_or_buf=outpath+sp+'_'+test_set+'_set.test.txt', sep='\t')


#with open('filename.pickle', 'wb') as handle:
#    pickle.dump(a, handle, protocol=pickle.HIGHEST_PROTOCOL)


# add species....
n_perms = 50000
cores = 6
all_test_obj = []
all_permutations = []
all_per_set = []
for test_set in tests:
    function_sets = psp.FunctionSets(function_set_file='../gene_sets/'+test_set+'_long_reformat.txt', min_set_size=5, annotation_obj=annotations)
    species = ['abelii', 'pygmaeus', 'gorilla', 'graueri', 'bonobo', 'troglodytes', 'schweinfurthii', 'ellioti', 'verus']
    for sp in species:
        candidates = psp.Variants(variant_file='../gene_set_test_files/'+sp+'_candidate_windows.txt')
        candidates.annotate_variants(annotation_obj=annotations)
        background = psp.Variants(variant_file='../gene_set_test_files/'+sp+'_background_windows.txt')
        background.annotate_variants(annotation_obj=annotations)
        test_obj = psp.TestObject(candidates,
                                  background,
                                  function_sets,
                                  n_cores=cores)
        all_test_obj.append(test_obj)
        permutations = psp.Permutation(test_obj, n_perms, cores)
        all_permutations.append(permutations)
        per_set = psp.SetPerPerm(permutations,
                                 function_sets,
                                 test_obj,
                                 cores)
        all_per_set.append(per_set)