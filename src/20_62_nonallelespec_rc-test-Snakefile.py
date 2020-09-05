
# Snakefile for running the Row-Column test
# for co-mutation and mutual exclusivity

from pathlib import Path
# import drmaa
import os

# permutation numbers
NUM_PERMS = 10 #10000
PERMUTATIONS = list(range(1, NUM_PERMS + 1))


data_dir = 'data/rc-test-nonallelespec/'
perm_gr_dir = '/n/scratch3/users/j/jc604/rc-test-nonallelespecific_permuted_graphs/'


wildcard_constraints:
	cancer='COAD|LUAD|MM|PAAD'


# run the pipeline
rule all:
	input:
		data_dir + 'output/exclusivity/COAD_exclusivity_results.rds',
		data_dir + 'output/exclusivity/LUAD_exclusivity_results.rds',
		data_dir + 'output/exclusivity/MM_exclusivity_results.rds',
		data_dir + 'output/exclusivity/PAAD_exclusivity_results.rds'




# create the bipartite graph from the table of samples and their mutated genes
rule make_bipartite_graph:
	input:
		in_file = data_dir + 'input/{cancer}_mutations.tsv'
	output:
		save_name = data_dir + 'intermediate/{cancer}_bgr.rds'
	script:
		'20_05_make-bipartite-graph.R'




perm_time_alloc_dict = {
	'COAD': '06:30:00',
	'LUAD': '24:00:00',
	  'MM': '01:00:00',
	'PAAD': '01:45:00',
	'SKCM': '48:00:00',
}
perm_mem_alloc_dict = {
	'COAD': '270',
	'LUAD': '250',
	  'MM': '240',
	'PAAD': '240',
	'SKCM': '300',
}

# use edge swapping to make the permutation matrices
rule permute_bipartite_graph:
	input:
		bipartite_gr = data_dir + 'intermediate/{cancer}_bgr.rds'
	output:
		save_name = perm_gr_dir + '{cancer}/{cancer}_perm{perm_num}.rds'
	params:
		time_alloc = lambda w: perm_time_alloc_dict[w.cancer],
		mem = lambda w: perm_mem_alloc_dict[w.cancer],
		partition = lambda w: 'medium' if w.cancer == 'LUAD' else 'short',
		Q = 20
	script:
		'20_06_permute-bipartite-graph.R'




results_df_time_dict = {
	'COAD': {'comutation': '00:35:00', 'exclusivity': '00:25:00'},
	'LUAD': {'comutation': '00:35:00', 'exclusivity': '00:25:00'},
	  'MM': {'comutation': '00:13:00', 'exclusivity': '00:08:00'},
	'PAAD': {'comutation': '00:10:00', 'exclusivity': '00:10:00'},
	'SKCM': {'comutation': '00:10:00', 'exclusivity': '00:10:00'},
}
results_df_mem_dict = {
	'COAD': {'comutation': 38000, 'exclusivity': 27000},
	'LUAD': {'comutation': 39000, 'exclusivity': 39000},
	  'MM': {'comutation': 20000, 'exclusivity': 10000},
	'PAAD': {'comutation': 30000, 'exclusivity': 20000},
	'SKCM': {'comutation': 20000, 'exclusivity': 20000},
}
rule make_results_df:
	input: 
		real_gr = data_dir + 'intermediate/{cancer}_bgr.rds'
	output:
		output_name = data_dir + 'intermediate/{cancer}_{which_test}_results_df.rds'
	params:
		which_test_short = lambda w: 'co' if w.which_test == 'comutation' else 'ex',
		time = lambda w: results_df_time_dict[w.cancer][w.which_test],
		mem = lambda w: results_df_mem_dict[w.cancer][w.which_test],
		min_times_mut = lambda w: 2 if w.which_test == 'comutation' else 3
	script:
		'20_61_nonallelespec_make-results-df.R'




process_permuted_graph_time_dict = {
	'COAD': {'comutation': '48:00:00', 'exclusivity': '200:00:00'},
	'LUAD': {'comutation': '200:00:00', 'exclusivity': '200:00:00'},
	  'MM': {'comutation': '48:00:00', 'exclusivity': '120:00:00'},
	'PAAD': {'comutation': '120:00:00', 'exclusivity': '120:00:00'},
	'SKCM': {'comutation': '48:00:00', 'exclusivity': '120:00:00'},
}
process_permuted_graph_mem_dict = {
	'COAD': {'comutation': 11000, 'exclusivity': 50000},
	'LUAD': {'comutation': 30000, 'exclusivity': 50000},
	  'MM': {'comutation': 11000, 'exclusivity': 35000},
	'PAAD': {'comutation': 35000, 'exclusivity': 35000},
	'SKCM': {'comutation': 11000, 'exclusivity': 35000},
}

rule process_permuted_graph:
	input:
		results_df = data_dir + 'intermediate/{cancer}_{which_test}_results_df.rds',
		perm_grs = expand(perm_gr_dir + '{{cancer}}/{{cancer}}_perm{perm_num}.rds',
		                  perm_num=PERMUTATIONS)
	output:
		save_name = data_dir + 'output/{which_test}/{cancer}_{which_test}_results.rds'
	params:
		n_cores = 20,
		which_test_short = lambda w: 'co' if w.which_test == 'comutation' else 'ex',
		time = lambda w: process_permuted_graph_time_dict[w.cancer][w.which_test],
		mem = lambda w: process_permuted_graph_mem_dict[w.cancer][w.which_test],
		partition = lambda w: 'medium' if w.which_test == 'comutation' else 'long'
	script:
		'20_08_process-permuted-graph.R'
