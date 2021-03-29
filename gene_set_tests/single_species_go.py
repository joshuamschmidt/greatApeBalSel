import pysetperm as psp
import numpy as np
import pandas as pd

n_perms = 50000
cores = 4
annotations = psp.AnnotationSet(annotation_file='data/genes.txt', range_modification=2000)
function_sets = psp.FunctionSets(function_set_file='data/kegg.txt', min_set_size=10, annotation_obj=annotations)
# specific inputs
e_candidates = psp.Variants(variant_file='data/eastern_candidates.txt')
e_candidates.annotate_variants(annotation_obj=annotations)
e_background = psp.Variants(variant_file='data/eastern_background.txt.gz')
e_background.annotate_variants(annotation_obj=annotations)
