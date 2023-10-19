# rule liger_lambda_and_k_sweep:
# 	input:
# 		scelist = dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNAseq-{subtype}-scelist-{condition}.rds',
# 	params:
# 		seed = '{seed}',
# 		Lambda = '{Lambda}',
# 		k = '{k}',
# 	output:
# 		metrics = resultoutput + 'LIGER/parameter-sweep/metrics/{subtype}/{condition}/seed-{seed}/{subtype}-{condition}-seed-{seed}-lambda-{Lambda}-k-{k}.tsv',
# 		factor_contributions = resultoutput + 'LIGER/parameter-sweep/factor-contributions/{subtype}/{condition}/seed-{seed}/k-{k}/{subtype}-{condition}-seed-{seed}-lambda-{Lambda}-k-{k}.tsv',
# 	resources:
# 		mem_mb = 8000
# 	threads: 1
# 	script:
# 		'signature-extraction/liger-lambda-and-k-sweep.R'

# rule combine_liger_parameter_sweep_metrics:
# 	input:
# 		tsv = expand(resultoutput + 'LIGER/parameter-sweep/metrics/{{subtype}}/{{condition}}/seed-{seed}/{{subtype}}-{{condition}}-seed-{seed}-lambda-{Lambda}-k-{k}.tsv', seed = seeds, Lambda = lambdalist, k = klist),
# 	params:
# 	output:
# 		tsv = resultoutput + 'LIGER/parameter-sweep/metrics/{subtype}/{condition}/{subtype}-{condition}-parameter-sweep-metrics.tsv',
# 	resources:
# 		mem_mb = 2000
# 	threads: 1
# 	script:
# 		'signature-extraction/combine-liger-parameter-sweep-metrics.R'

# rule visualize_liger_parameter_sweep_metrics:
# 	input:
# 		tsv = resultoutput + 'LIGER/parameter-sweep/metrics/{subtype}/{condition}/{subtype}-{condition}-parameter-sweep-metrics.tsv',
# 	params:
# 	output:
# 		plot_k = figureoutput + 'LIGER/parameter-sweep/metrics/{subtype}/{condition}/{subtype}-{condition}-k-sweep.png',
# 		plot_lambda = figureoutput + 'LIGER/parameter-sweep/metrics/{subtype}/{condition}/{subtype}-{condition}-lambda-sweep.png',
# 	resources:
# 		mem_mb = 1000
# 	threads: 1
# 	script:
# 		'signature-extraction/visualize-liger-parameter-sweep-metrics.R'

rule liger_variable_gene_selection:
	input:
		sces_dis = dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNASeq-{subtype}-scelist-discovery.rds',
		sces_val = dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNASeq-{subtype}-scelist-validation.rds',
	params:
		exprs_var_thres = config['liger']['exprssion_variation_threshold'],
	output:
		var_genes_dis = resultoutput + 'LIGER/signature-extraction/variable-gene-list/{subtype}/{subtype}-variable-genes-discovery.tsv',
		var_genes_val = resultoutput + 'LIGER/signature-extraction/variable-gene-list/{subtype}/{subtype}-variable-genes-validation.tsv',
		var_genes_common = resultoutput + 'LIGER/signature-extraction/variable-gene-list/{subtype}/{subtype}-variable-genes-common.tsv',
	resources:
		mem_mb = 10000
	threads: 4
	script:
		'signature-extraction/liger-variable-gene-selection.R'

rule liger_signature_extraction:
	input:
		scelist = dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNASeq-{subtype}-scelist-{condition}.rds',
		var_genes_common = resultoutput + 'LIGER/signature-extraction/variable-gene-list/{subtype}/{subtype}-variable-genes-common.tsv',
	params:
		k = lambda wildcards: get_liger_param(parameterlist, wildcards.subtype, 'k'),
		Lambda = lambda wildcards: get_liger_param(parameterlist, wildcards.subtype, 'lambda'),
		seed = lambda wildcards: get_liger_param(parameterlist, wildcards.subtype, 'seed'),
		# liger = resultoutput + 'LIGER/signature-extraction/LIGER-object/{subtype}/{subtype}-liger-{condition}.rds',
	output:
		liger = resultoutput + 'LIGER/signature-extraction/LIGER-object/{subtype}/{subtype}-liger-{condition}.rds',
	resources:
		mem_mb = 50000
	threads: 16
	script:
		'signature-extraction/liger-signature-extraction.R'

rule visualize_liger_output:
	input:
		liger = resultoutput + 'LIGER/signature-extraction/LIGER-object/{subtype}/{subtype}-liger-{condition}.rds',
	params:
		gene_loadings_plot_dir = figureoutput + 'LIGER/signature-extraction/{subtype}/{condition}/gene-loadings/',
	output:
		signature_loading_plot = figureoutput + 'LIGER/signature-extraction/{subtype}/{condition}/{subtype}-{condition}-signature-loading.png',
		gene_loading_plot = figureoutput + 'LIGER/signature-extraction/{subtype}/{condition}/{subtype}-{condition}-gene-loading.png',
		dataset_cluster_plot = figureoutput + 'LIGER/signature-extraction/{subtype}/{condition}/{subtype}-{condition}-dataset-cluster.png',
	resources:
		mem_mb = 8000
	threads: 4
	script:
		'signature-extraction/visualize-liger-output.R'

rule extract_liger_loading_matrices:
	input:
		liger_dis = resultoutput + 'LIGER/signature-extraction/LIGER-object/{subtype}/{subtype}-liger-discovery.rds',
		liger_val = resultoutput + 'LIGER/signature-extraction/LIGER-object/{subtype}/{subtype}-liger-validation.rds',
		sce_dis = dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNASeq-{subtype}-sce-discovery.rds',
		sce_val = dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNASeq-{subtype}-sce-validation.rds',
	params:
	output:
		gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading.tsv',
		sig_loading_mtx_dis = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-discovery.tsv',
		sig_loading_mtx_val = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-validation.tsv',
	resources:
		mem_mb = 8000
	threads: 1
	script:
		'signature-extraction/extract-liger-loading-matrices.R'

