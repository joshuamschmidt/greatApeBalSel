import pysetperm as psp
import numpy as np
import pandas as pd

n_perms = 50000
cores = 4
annotations = psp.AnnotationSet(annotation_file='../gene_annotations/hg19_nonRedundant_GOWINDA_noORS.bed', range_modification=2000)
function_sets = psp.FunctionSets(function_set_file='../gene_sets/kegg_long.txt', min_set_size=10, annotation_obj=annotations)
# specific inputs
# eastern chimps
e_candidates = psp.Variants(variant_file='../gene_set_test_files/schweinfurthii_candidate_windows.txt')
e_candidates.annotate_variants(annotation_obj=annotations)
e_background = psp.Variants(variant_file='../gene_set_test_files/schweinfurthii_background_windows.txt')
e_background.annotate_variants(annotation_obj=annotations)
e_test_obj = psp.TestObject(e_candidates,
                            e_background,
                            function_sets,
                            n_cores=cores)
e_permutations = psp.Permutation(e_test_obj, n_perms, cores)
e_per_set = psp.SetPerPerm(e_permutations,
                           function_sets,
                           e_test_obj,
                           cores)
e_results = make_results_table(e_test_obj, function_sets, e_per_set)
