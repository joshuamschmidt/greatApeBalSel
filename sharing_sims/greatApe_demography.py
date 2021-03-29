import math
import msprime
import itertools
import allel
import numpy as np


def great_ape_demography(num_replicates, sim_locus_length):
    n_ms: int = 10000
    generation_time: int = 20
    m_scale = 4 * n_ms
    # First we set out the values for pop sizes.
    n_present_human = 10000
    n_present_central = 30000
    n_present_eastern = 12000
    n_present_nigeria = 9000
    n_present_western = 5000
    n_present_bonobo = 5000
    n_present_wlg = 15000
    n_present_elg = 2000
    n_present_bo = 8000
    n_present_so = 25000
    n_present_macaque = 35000
    # growth rates from present to some past date....
    # size change at 100 kya
    t_100kya: int = 100e3 / generation_time
    n_100kya_western = 30000
    n_100kya_nigeria = 15000
    # size change at 150 kya
    t_150kya: int = 150e3 / generation_time
    n_150kya_bonobo = 20000
    # size change at 200 kya
    n_200kya_central = 40000
    n_200kya_eastern = 30000
    n_200kya_elg = 15000
    # pop splits at 200 kya central/eastern; wlg/elg
    t_200kya: int = 200e3 / generation_time
    # and sizes of joint lineages
    n_central_eastern_split = 30000
    n_wlg_elg_split = 15000
    # size changes and pop splits at 300 kya; western nigeria split.
    t_300kya: int = 300e3 / generation_time
    n_300kya_western = 5000
    n_300kya_nigeria = 10000
    n_300kya_so = 40000
    n_nigeria_western_split = 10000
    # at 500 kya: orangs split; all chimps split. adjust pop sizes and growth rates.
    t_500kya: int = 500e3 / generation_time
    n_500kya_bo = 20000
    n_500kya_so = 20000
    n_500kya_central_eastern = 20000
    n_500kya_nigeria_western = 15000
    n_500kya_orangs = 20000
    n_500kya_all_chimps = 20000
    # size change at 600 kya
    t_600kya: int = 600e3 / generation_time
    n_600kya_bonobo = 5000
    # size change at 1mya
    t_1Mya: int = 1000e3 / generation_time
    n_1mya_bonobo = 10000
    n_1mya_all_chimps = 10000
    n_1mya_wlg_elg = 20000
    n_all_pan_split = 10000
    # size change at 2.5mya
    t_2_5Mya: int = 2500e3 / generation_time
    n_2_5mya_wlg_elg = 12000
    # pan and homo join at 5Mya
    t_5Mya: int = 5000e3 / generation_time
    n_5Mya_all_pan = 23469
    n_pan_human_split = 23469
    # size change at 7.5mya
    t_7_5Mya: int = 7500e3 / generation_time
    n_7_5mya_all_pan = 40000
    n_7_5mya_wlg_elg = 30000
    n_african_ape_split = 40000
    # size change at 12.5mya
    t_12_5Mya: int = 12500e3 / generation_time
    n_12_5mya_african_ape = 60000
    n_12_5mya_orangs = 20000
    n_great_ape_split = 125000
    # size change at 25mya
    t_25Mya: int = 25000e3 / generation_time
    n_25mya_great_ape = 125000
    n_25mya_macaque = 35000
    n_ape_macaque_split = 125000
    # didnt notice the difference in n_present_western as unscaled this is still just 0.081.
    # what re the demes?
    # 0 human
    # 1 central chimp
    # 2 eastern chimp
    # 3 nigeria chimp
    # 4 western chimp
    # 5 bonobo
    # 6 wl gorilla
    # 7 el gorilla
    # 8 b orang
    # 9 s orang
    # 10 macaque
    # for growth rates use my R function gRatePercentage
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=1, initial_size=n_present_human),
        msprime.PopulationConfiguration(
            sample_size=8, initial_size=n_present_central, growth_rate=-0.00002876821),
        msprime.PopulationConfiguration(
            sample_size=12, initial_size=n_present_eastern, growth_rate=-0.00009162907),
        msprime.PopulationConfiguration(
            sample_size=20, initial_size=n_present_nigeria, growth_rate=-0.0001021651),
        msprime.PopulationConfiguration(
            sample_size=8, initial_size=n_present_western, growth_rate=-0.0003583519),
        msprime.PopulationConfiguration(
            sample_size=26, initial_size=n_present_bonobo, growth_rate=-0.0001848392),
        msprime.PopulationConfiguration(
            sample_size=46, initial_size=n_present_wlg),
        msprime.PopulationConfiguration(
            sample_size=6, initial_size=n_present_elg, growth_rate=-0.0002014903),
        msprime.PopulationConfiguration(
            sample_size=10, initial_size=n_present_bo, growth_rate=-0.00003665163),
        msprime.PopulationConfiguration(
            sample_size=10, initial_size=n_present_so, growth_rate=-0.00003133358),
        msprime.PopulationConfiguration(
            sample_size=1, initial_size=n_present_macaque),
    ]
    migration_matrix = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    ]
    demographic_events = [
        # define all demographic events into the past....
        # size change at 100 kya for nigeria and western. also adjust growth rates back to 300kya
        msprime.PopulationParametersChange(
            time=t_100kya, initial_size=n_100kya_nigeria, growth_rate=0.00004054651, population_id=3),
        msprime.PopulationParametersChange(
            time=t_100kya, initial_size=n_100kya_western, growth_rate=0.0001791759, population_id=4),
        # growth rate change for bonobo at 150kya back to 600kya, 20k to 5k
        msprime.PopulationParametersChange(
            time=t_150kya, initial_size=n_150kya_bonobo, growth_rate=0.00006161308, population_id=5),
        # at 200kya - central and eastern split. wlg and elg split. set new Ne and growth rates too.
        msprime.PopulationParametersChange(
            time=t_200kya, initial_size=n_200kya_central, growth_rate=0.0, population_id=1),
        msprime.PopulationParametersChange(
            time=t_200kya, initial_size=n_200kya_eastern, growth_rate=0.0, population_id=2),
        msprime.PopulationParametersChange(
            time=t_200kya, initial_size=n_200kya_elg, growth_rate=0.0, population_id=7),
        msprime.MassMigration(
            time=t_200kya, source=2, destination=1, proportion=1.0),
        msprime.MassMigration(
            time=t_200kya, source=7, destination=6, proportion=1.0),
        msprime.PopulationParametersChange(
            time=t_200kya, initial_size=n_central_eastern_split, growth_rate=0.00002703101, population_id=1),
        msprime.PopulationParametersChange(
            time=t_200kya, initial_size=n_wlg_elg_split, growth_rate=-0.000007192052, population_id=6),
        # at 300kya - nigeria and western split. set new Ne and growth rates to 500kya. S orang rate change to 500Kya
        msprime.PopulationParametersChange(
            time=t_300kya, initial_size=n_300kya_nigeria, growth_rate=0.0, population_id=3),
        msprime.PopulationParametersChange(
            time=t_300kya, initial_size=n_300kya_western, growth_rate=0.0, population_id=4),
        msprime.MassMigration(
            time=t_300kya, source=4, destination=3, proportion=1.0),
        msprime.PopulationParametersChange(
            time=t_300kya, initial_size=n_nigeria_western_split, growth_rate=-0.00004054651, population_id=3),
        msprime.PopulationParametersChange(
            time=t_300kya, initial_size=n_300kya_so, growth_rate=0.00006931472, population_id=9),
        # at 500 kya: orangs split; all chimps split. adjust pop sizes and growth rates.
        msprime.PopulationParametersChange(
            time=t_500kya, initial_size=n_500kya_central_eastern, growth_rate=0.0, population_id=1),
        msprime.PopulationParametersChange(
            time=t_500kya, initial_size=n_500kya_nigeria_western, growth_rate=0.0, population_id=3),
        msprime.PopulationParametersChange(
            time=t_500kya, initial_size=n_500kya_bo, growth_rate=0.0, population_id=8),
        msprime.PopulationParametersChange(
            time=t_500kya, initial_size=n_500kya_so, growth_rate=0.0, population_id=9),
        msprime.MassMigration(
            time=t_500kya, source=3, destination=1, proportion=1.0),
        msprime.MassMigration(
            time=t_500kya, source=9, destination=8, proportion=1.0),
        msprime.PopulationParametersChange(
            time=t_500kya, initial_size=n_500kya_all_chimps, growth_rate=0.00002772589, population_id=1),
        msprime.PopulationParametersChange(
            time=t_500kya, initial_size=n_500kya_orangs, growth_rate=0, population_id=8),
        # at 600 kya bonobo growth rate change back to 1MYA
        msprime.PopulationParametersChange(
            time=t_600kya, initial_size=n_600kya_bonobo, growth_rate=-0.00003465736, population_id=5),
        # at 1MYA: chimps and bonobo join. chnage gorilla growth rate.
        msprime.PopulationParametersChange(
            time=t_1Mya, initial_size=n_1mya_all_chimps, growth_rate=0, population_id=1),
        msprime.PopulationParametersChange(
            time=t_1Mya, initial_size=n_1mya_bonobo, growth_rate=0, population_id=5),
        msprime.MassMigration(
            time=t_1Mya, source=5, destination=1, proportion=1.0),
        msprime.PopulationParametersChange(
            time=t_1Mya, initial_size=n_all_pan_split, growth_rate=-0.000004265521, population_id=1),
        msprime.PopulationParametersChange(
            time=t_1Mya, initial_size=n_1mya_wlg_elg, growth_rate=0.000006811008, population_id=6),
        # at 2.5MYA: change gorilla growth rate back to 7.5 MYA.
        msprime.PopulationParametersChange(
            time=t_2_5Mya, initial_size=n_2_5mya_wlg_elg, growth_rate=-0.000003665163, population_id=6),
        # at 5MYA: pan and homo split. keep Pan back in time i.e. 0 joins 1.
        msprime.PopulationParametersChange(
            time=t_5Mya, initial_size=n_5Mya_all_pan, growth_rate=0, population_id=1),
        msprime.MassMigration(
            time=t_5Mya, source=0, destination=1, proportion=1.0),
        msprime.PopulationParametersChange(
            time=t_5Mya, initial_size=n_pan_human_split, growth_rate=-0.000004265521, population_id=1),
        # at 7.5 MYya. size changes. pan/homo and gorilla split
        msprime.PopulationParametersChange(
            time=t_7_5Mya, initial_size=n_7_5mya_all_pan, growth_rate=0, population_id=1),
        msprime.PopulationParametersChange(
            time=t_7_5Mya, initial_size=n_7_5mya_wlg_elg, growth_rate=0, population_id=6),
        msprime.MassMigration(
            time=t_7_5Mya, source=6, destination=1, proportion=1.0),
        msprime.PopulationParametersChange(
            time=t_7_5Mya, initial_size=n_african_ape_split, growth_rate=-0.00000162186, population_id=1),
        # 12.5 Mya African apes and Orangs split.
        msprime.PopulationParametersChange(
            time=t_12_5Mya, initial_size=n_12_5mya_african_ape, growth_rate=0, population_id=1),
        msprime.PopulationParametersChange(
            time=t_12_5Mya, initial_size=n_12_5mya_orangs, growth_rate=0, population_id=8),
        msprime.MassMigration(
            time=t_12_5Mya, source=8, destination=1, proportion=1.0),
        msprime.PopulationParametersChange(
            time=t_12_5Mya, initial_size=n_great_ape_split, growth_rate=0, population_id=1),
        # 25 Mya Great Apes Macaqsa split
        msprime.PopulationParametersChange(
            time=t_25Mya, initial_size=n_25mya_great_ape, growth_rate=0, population_id=1),
        msprime.PopulationParametersChange(
            time=t_25Mya, initial_size=n_25mya_macaque, growth_rate=0, population_id=10),
        msprime.MassMigration(
            time=t_25Mya, source=10, destination=1, proportion=1.0),
        msprime.PopulationParametersChange(
            time=t_25Mya, initial_size=n_ape_macaque_split, growth_rate=0, population_id=1)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    # dd = msprime.DemographyDebugger(
    #     population_configurations=population_configurations,
    #     demographic_events=demographic_events,
    #     migration_matrix=None)
    # dd.print_history()
    tree_seq: object = msprime.simulate(
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        length=sim_locus_length,
        recombination_rate=1e-8,
        mutation_rate=1e-9,
        num_replicates=num_replicates
    )
    return tree_seq


sample_config = [1, 8, 12, 20, 8, 26, 46, 6, 10, 10, 1]
locus_length = 30e3
window_length = 3000
window_slide = 1500

sliding_windows_starts = np.arange(1, locus_length, window_slide)


def tree():
    tree_seq = great_ape_demography(1, 30000)
    for a, tree_sequence in enumerate(tree_seq):
        if a == 0:
            first_tree = tree_sequence.first()
            print(first_tree.draw(format="unicode"))


rep_ts = great_ape_demography(1, sim_locus_length=30000)

for j, ts in enumerate(rep_ts):
    print(j)
    # get the index for the genomes sampled from each population
    # note we can create a list with list comprehension (the pythonic way)
    pop_list = [pop.id for pop in ts.populations()]
    pop_indices = [ts.samples(population=pop) for pop in pop_list]
    # now lets get the variant information
    V = np.zeros((ts.get_num_mutations(), ts.get_sample_size()),
                 dtype=np.int8)
    for variant in ts.variants():
        V[variant.index] = variant.genotypes
    # create a haplotype array in allele from numpy arrays of 0s/1s
    # as the dataset is reasonably small there is no need to use the chunking functionality of allel
    haplotypes = allel.HaplotypeArray(V)
    # create a list of allele counts for each population
    allele_counts = [haplotypes.count_alleles(max_allele=None, subpop=pop_idx) for pop_idx in pop_indices]
    ## calculate the hudson fsts
    all_fsts = np.zeros((len(pop_list), len(pop_list)),
                        dtype=np.float64)
    for pair in itertools.combinations(pop_list, 2):
        num, den = allel.hudson_fst(allele_counts[pair[0]], allele_counts[pair[1]])
        all_fsts[pair[1], pair[0]] = np.sum(num) / np.sum(den)
    print(all_fsts)

for j, ts in enumerate(rep_ts):
    print(j)
    for variant in ts.variants():
        print(
            variant.site.id, variant.site.position,
            variant.alleles, variant.genotypes, sep="\t")

# get all the allele counts.
sample_config = [1, 8, 12, 20, 8, 26, 46, 6, 10, 10, 1]
genotypes = ts.genotype_matrix()
pops = [pop.id for pop in ts.populations()]
pop_indices = [ts.samples(population=pop) for pop in pops]
# create a haplotype array in allele from numpy arrays of 0s/1s
# as the dataset is reasonably small there is no need to use the chunking functionality of allele
haplotypes = allel.HaplotypeArray(genotypes)
# create a list of allele counts for each population
allele_counts = [haplotypes.count_alleles(max_allele=None, subpop=pop_idx) for pop_idx in pop_indices]
# convert to frequencies.
allele_frequencies = []
for j, nchr in enumerate(sample_config):
    allele_frequencies.append(allele_counts[j] / nchr)

allele_counts[2] / 12
# need to keep track of sites that are duplicated.....
all_sites = []
for j, site in enumerate(ts.sites()):
    physpos = int(round(site.position))
    if physpos == 0:
        physpos = 1
    while physpos in all_sites:
        physpos = physpos + 1
    all_sites.append(physpos)

dsite = []
for j, site in enumerate(ts.sites()):
    if j == 1:
        dsite = site


def calc_ncd(pop_idxs, outgroup_idx, target_frequency):

    return ncd
