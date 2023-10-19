rule liger_signature_validation:
	input:
		gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading.tsv',
	params:
		corr_method = config['signatures']['corr_method'],
		corr_thres = config['signatures']['corr_thres'],
		dist_method = config['signatures']['dist_method'],
		dist_thres = config['signatures']['dist_thres'],
		minkowski_p = config['signatures']['minkowski_p']
	output:
		gene_loading_corr = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-corr.tsv',
		gene_loading_dist = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-dist.tsv',
		validated_sig_df = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-validated-signatures.tsv',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'signature-validation/liger-signature-validation.R'

rule visualize_liger_signature_validation:
	input:
		gene_loading_corr = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-corr.tsv',
		gene_loading_dist = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-dist.tsv',
	params:
		dist_method = config['signatures']['dist_method'],
	output:
		gene_loading_corr_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-corr.png',
		gene_loading_dist_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-dist.png',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'signature-validation/visualize-liger-signature-validation.R'

rule liger_signature_validation_alternative_methods:
	input:
		gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading.tsv',
	params:
		corr_method = config['signatures']['corr_method'],
		corr_thres = config['signatures']['corr_thres'],
		dist_method = config['signatures']['dist_method'],
		dist_thres = config['signatures']['dist_thres'],
		minkowski_p = config['signatures']['minkowski_p']
	output:
		#gene_loading_corr = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-corr.tsv',
		#gene_loading_dist = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-dist.tsv',
		validated_sig_df = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-validated-signatures-hungarian-method.tsv',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'signature-validation/liger-signature-validation-alternative-methods.R'

rule visualize_liger_signature_validation_alternative_methods:
	input:
		gene_loading_corr = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-corr.tsv',
		gene_loading_dist = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-dist.tsv',
		validated_sig_df = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-validated-signatures-hungarian-method.tsv',
	params:
		dist_method = config['signatures']['dist_method'],
	output:
		gene_loading_corr_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-corr-hungarian-method.png',
		gene_loading_dist_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-dist-hungarian-method.png',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'signature-validation/visualize-liger-signature-validation-alternative-methods.R'

rule liger_signature_validation_simulation_test:
	input:
		gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading.tsv',
	params:
		corr_method = config['signatures']['corr_method'],
		#corr_thres = config['signatures']['corr_thres'],
		dist_method = config['signatures']['dist_method'],
		#dist_thres = config['signatures']['dist_thres'],
		minkowski_p = config['signatures']['minkowski_p'],
		num_sim = config['signatures']['num_sim'],
	output:
		simulated_dist_list = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-simulated-distance-list.tsv',
		simulated_corr_list = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-simulated-correlation-list.tsv',
	resources:
		mem_mb = 10000
	threads: 2
	script:
		'signature-validation/liger-signature-validation-simulation-test.R'

rule visualize_liger_signature_validation_simulation_test:
	input:
		simulated_dist_list = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-simulated-distance-list.tsv',
		simulated_corr_list = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-simulated-correlation-list.tsv',
		gene_loading_dist = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-dist.tsv',
		gene_loading_corr = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-corr.tsv',

	params:
		#dist_method = config['signatures']['dist_method'],
	output:
		gene_loading_sim_dist_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-simulated-distance.pdf',
		gene_loading_sim_corr_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-simulated-correlation.pdf',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'signature-validation/visualize-liger-signature-validation-simulation-test.R'

rule extract_validated_liger_signatures:
	input:
		sig_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-discovery.tsv',
		gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading.tsv',
		validated_sig_df = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-validated-signatures.tsv',
	params:
		validation_metric = config['signatures']['validation_metric_to_use']
	output:
		validated_sig_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-validated.tsv',
		validated_gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-validated.tsv',
	resources: 
		mem_mb = 1000
	threads: 1
	script:
		'signature-validation/extract-validated-liger-signatures.R'

rule analyze_validated_liger_signatures:
	input:
		sig_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-validated.tsv',
		gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-validated.tsv',
	params:
	output:
		validated_sig_loading_corr = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-collapse/{subtype}-signature-loading-correlation-validated.tsv',
		validated_gene_loading_corr = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-collapse/{subtype}-gene-loading-correlation-validated.tsv',
	resources: 
		mem_mb = 1000
	threads: 1
	script:
		'signature-validation/analyze-validated-liger-signatures.R'

rule visualize_validated_liger_signatures_analysis:
	input:
		validated_sig_loading_corr = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-collapse/{subtype}-signature-loading-correlation-validated.tsv',
		validated_gene_loading_corr = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-collapse/{subtype}-gene-loading-correlation-validated.tsv',
	params:
	output:
		validated_gene_loading_corr_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-collapse/{subtype}-gene-loading-correlation-validated.png',
		validated_sig_loading_corr_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-collapse/{subtype}-signature-loading-correlation-validated.png',
	resources: 
		mem_mb = 1000
	threads: 1
	script:
		'signature-validation/visualize-validated-liger-signatures-analysis.R'

rule generate_liger_signature_collapse_guide:
	input:
		#sig_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-validated.tsv',
		#gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-validated.tsv',
		validated_gene_loading_corr = lambda wildcards: expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-collapse/{subtype}-gene-loading-correlation-validated.tsv', subtype = ct_extract),
	params:
		collapse_corr_threshold = config['signatures']['collapse_corr_threshold'],
	output:
		#validated_sig_loading_corr = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-collapse/{subtype}-signature-loading-correlation-validated.tsv',
		#validated_gene_loading_corr = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-collapse/{subtype}-gene-loading-correlation-validated.tsv',
		signature_collapse_guide = resultoutput + 'LIGER/signature-analysis/signature-collapse-guide.tsv',
	resources: 
		mem_mb = 1000
	threads: 1
	script:
		'signature-validation/generate-liger-signature-collapse-guide.R'

rule liger_signature_collapse:
	input:
		gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-validated.tsv',
		sig_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-validated.tsv',
		signature_collapse_guide = resultoutput + 'LIGER/signature-analysis/signature-collapse-guide.tsv',
	params:
		
	output:
		collapsed_sig_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-collapsed.tsv',
		collapsed_gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-collapsed.tsv',
		sig_collapse_namemap = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-collapse-name-mapping.tsv',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'signature-validation/liger-signature-collapse.R'

rule liger_score_validation_with_collapsed_signatures:
	input:
		gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-collapsed.tsv',
		scelist_val = dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNASeq-{subtype}-scelist-validation.rds',
	params:
		score_each_validation_cohort = config['signatures']['score_each_validation_cohort'],
		assay_used_for_scoring = config['signatures']['assay_used_for_scoring'],
		cohort_correct_weight = config['signatures']['cohort_correct_weight'],
	output:
		scored_sig_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-collapsed-scored-validation.tsv',
		scored_gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-collapsed-scored-validation.tsv',
		#sig_score_namemap = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-score-name-mapping.tsv',
	resources:
		mem_mb = 20000
	threads: 10
	script:
		'signature-validation/liger-score-validation-with-collapsed-signatures.R'

