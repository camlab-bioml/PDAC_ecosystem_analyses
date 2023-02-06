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
		'signature-analysis/liger-signature-validation.R'

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
		'signature-analysis/visualize-liger-signature-validation.R'

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
		'signature-analysis/extract-validated-liger-signatures.R'

rule liger_signature_collapse:
	input:
		gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-validated.tsv',
		sig_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-validated.tsv',
		signature_collapse_guide = config['signatures']['signature_collapse_guide'],
	params:
		
	output:
		collapsed_sig_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-collapsed.tsv',
		collapsed_gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-collapsed.tsv',
		sig_collapse_namemap = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-collapse-name-mapping.tsv',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'signature-analysis/liger-signature-collapse.R'

rule liger_signature_top_gene_loading_analysis_disval:
	input:
		gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading.tsv',
	params:
		num_top_genes = config['signatures']['num_top_genes'],
	output:
		sig_top_gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-loading.tsv',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'signature-analysis/liger-signature-top-gene-loading-analysis-disval.R'

rule visualize_liger_signature_top_gene_loading_analysis_disval:
	input:
		sig_top_gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-loading.tsv',
	params:
		
	output:
		sig_top_gene_loading_heatmap = figureoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-loading.png',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'signature-analysis/visualize-liger-signature-top-gene-loading-analysis-disval.R'

rule liger_signature_top_gene_loading_analysis:
	input:
		gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-{condition}.tsv',
	params:
		num_top_genes = config['signatures']['num_top_genes'],
	output:
		sig_top_gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-loading-{condition}.tsv',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'signature-analysis/liger-signature-top-gene-loading-analysis.R'

rule visualize_liger_signature_top_gene_loading_analysis:
	input:
		sig_top_gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-loading-{condition}.tsv',
	params:
		
	output:
		sig_top_gene_loading_heatmap = figureoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-loading-{condition}.png',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'signature-analysis/visualize-liger-signature-top-gene-loading-analysis.R'

rule liger_signature_known_markers_analysis_disval:
	input:
		gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading.tsv',
		marker_list = config['signatures']['known_signature_markers_dir'] + '{subtype}/{marker_list}' + config['signatures']['known_signature_markers_file_pattern'],
	params:
	output:
		sig_known_gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/known-markers-loading/{subtype}-{marker_list}-gene-loading.tsv',
	resources: 
		mem_mb = 1000
	threads: 1
	script:
		'signature-analysis/liger-signature-known-markers-analysis-disval.R'

rule visualize_liger_signature_known_markers_analysis_disval:
	input:
		sig_known_gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/known-markers-loading/{subtype}-{marker_list}-gene-loading.tsv',
	params:
	output:
		sig_known_gene_loading_plot_full = figureoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/known-markers-loading/{subtype}-{marker_list}-gene-loading-full.png',
		sig_known_gene_loading_plot_cleaned = figureoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/known-markers-loading/{subtype}-{marker_list}-gene-loading-cleaned.png',
	resources: 
		mem_mb = 1000
	threads: 1
	script:
		'signature-analysis/visualize-liger-signature-known-markers-analysis-disval.R'

rule liger_signature_known_markers_analysis:
	input:
		gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-{condition}.tsv',
		marker_list = config['signatures']['known_signature_markers_dir'] + '{subtype}/{marker_list}' + config['signatures']['known_signature_markers_file_pattern'],
	params:
		num_top_genes_for_overlap_test = config['signatures']['num_top_genes_for_overlap_test'],
	output:
		sig_known_gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/known-markers-loading/{subtype}-{marker_list}-gene-loading-{condition}.tsv',
		sig_gom = resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/known-markers-overlap/{subtype}-{marker_list}-geneoverlap-matrix-{condition}.rds',
	resources: 
		mem_mb = 1000
	threads: 1
	script:
		'signature-analysis/liger-signature-known-markers-analysis.R'

rule visualize_liger_signature_known_markers_analysis:
	input:
		sig_known_gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/known-markers-loading/{subtype}-{marker_list}-gene-loading-{condition}.tsv',
		sig_gom = resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/known-markers-overlap/{subtype}-{marker_list}-geneoverlap-matrix-{condition}.rds',
	params:
		pvalue_cutoff_for_overlap_test = config['signatures']['pvalue_cutoff_for_overlap_test'],
	output:
		sig_known_gene_loading_plot_full = figureoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/known-markers-loading/{subtype}-{marker_list}-gene-loading-full-{condition}.png',
		sig_known_gene_loading_plot_cleaned = figureoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/known-markers-loading/{subtype}-{marker_list}-gene-loading-cleaned-{condition}.png',
		sig_gom_plot_oddsratio = figureoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/known-markers-overlap/{subtype}-{marker_list}-geneoverlap-oddsratio-{condition}.png',
		sig_gom_plot_jaccardindex = figureoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/known-markers-overlap/{subtype}-{marker_list}-geneoverlap-jaccardindex-{condition}.png',
	resources: 
		mem_mb = 1000
	threads: 1
	script:
		'signature-analysis/visualize-liger-signature-known-markers-analysis.R'

rule liger_signature_enrichment_analysis_extract_geneuniverse:
	input:
		scelist_dis = dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNAseq-{subtype}-scelist-discovery.rds',
	params:
		
	output:
		geneuniverse = resultoutput + 'LIGER/signature-analysis/{subtype}/enrichment-analysis/{subtype}-signature-geneuniverse-{condition}.rds',
	resources:
		mem_mb = 10000
	threads: 2
	script:
		'signature-analysis/liger-signature-enrichment-analysis-extract-geneuniverse.R'

rule liger_signature_enrichment_analysis:
	input:
		gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading.tsv',
		gene_loading_mtx_validated = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-validated.tsv',
		gene_loading_mtx_collapsed = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-collapsed.tsv',
		geneuniverse = resultoutput + 'LIGER/signature-analysis/{subtype}/enrichment-analysis/{subtype}-signature-geneuniverse-{condition}.rds',
	params:
		signature = '{k}',
		gsea_maxGSSize = config['signatures']['gsea_maxGSSize'],
	output:
		overrepresentation_go = resultoutput + 'LIGER/signature-analysis/{subtype}/enrichment-analysis/overrepresentation-analysis/{condition}/{subtype}-signature-{k}-overrepresentation-GO.rds',
		overrepresentation_kegg = resultoutput + 'LIGER/signature-analysis/{subtype}/enrichment-analysis/overrepresentation-analysis/{condition}/{subtype}-signature-{k}-overrepresentation-KEGG.rds',
		gsea_go = resultoutput + 'LIGER/signature-analysis/{subtype}/enrichment-analysis/GSEA/{condition}/{subtype}-signature-{k}-GSEA-GO.rds',
	resources:
		mem_mb = 8000
	threads: 4
	script:
		'signature-analysis/liger-signature-enrichment-analysis.R'

rule visualize_liger_signature_enrichment_analysis:
	input:
		overrepresentation_go = resultoutput + 'LIGER/signature-analysis/{subtype}/enrichment-analysis/overrepresentation-analysis/{condition}/{subtype}-signature-{k}-overrepresentation-GO.rds',
		overrepresentation_kegg = resultoutput + 'LIGER/signature-analysis/{subtype}/enrichment-analysis/overrepresentation-analysis/{condition}/{subtype}-signature-{k}-overrepresentation-KEGG.rds',
		gsea_go = resultoutput + 'LIGER/signature-analysis/{subtype}/enrichment-analysis/GSEA/{condition}/{subtype}-signature-{k}-GSEA-GO.rds',
		gene_loading_mtx_validated = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-validated.tsv',
		gene_loading_mtx_collapsed = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-collapsed.tsv',
	params:
		signature = '{k}',
	output:
		overrepresentation_go_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/enrichment-analysis/overrepresentation-analysis/{condition}/{subtype}-signature-{k}-overrepresentation-GO.png',
		overrepresentation_kegg_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/enrichment-analysis/overrepresentation-analysis/{condition}/{subtype}-signature-{k}-overrepresentation-KEGG.png',
		gsea_go_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/enrichment-analysis/GSEA/{condition}/{subtype}-signature-{k}-GSEA-GO.png',
	resources: 
		mem_mb = 1000
	threads: 1
	script:
		'signature-analysis/visualize-liger-signature-enrichment-analysis.R'

rule liger_signature_loading_analysis:
	input:
		sig_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-{condition}.tsv',
	params:
		sig_activation_threshold_quantile = config['signatures']['sig_activation_threshold_quantile']
	output:
		sig_loading_mtx_long = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-long-form-{condition}.tsv',
		sig_loading_uniqueness = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-uniqueness-{condition}.tsv',
		sig_loading_top_two = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-top-two-{condition}.tsv',
		sig_top_two_count = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-top-two-count-{condition}.tsv',
		sig_top_freq = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-top-frequency-{condition}.tsv',
		sig_loading_patient_summary = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-patient-summary-{condition}.tsv',
		sig_activation_freq = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-activation-frequency-{condition}.tsv',
	resources:
		mem_mb = 6000
	threads: 1
	script:
		'signature-analysis/liger-signature-loading-analysis.R'

rule visualize_liger_signature_loading_analysis:
	input:
		sig_loading_mtx_long = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-long-form-{condition}.tsv',
		sig_loading_uniqueness = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-uniqueness-{condition}.tsv',
		sig_loading_top_two = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-top-two-{condition}.tsv',
		sig_top_two_count = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-top-two-count-{condition}.tsv',
		sig_top_freq = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-top-frequency-{condition}.tsv',
	params:
		
	output:
		sig_loading_qqnorm_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-qqnorm-{condition}.png',
		sig_loading_auc_density_grid_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-auc-density-grid-{condition}.png',
		sig_loading_auc_density_grouped_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-auc-density-grouped-{condition}.png',
		sig_loading_auc_box_grouped_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-auc-box-grouped-{condition}.png',
		sig_loading_auc_box_summary_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-auc-box-summary-{condition}.png',
		sig_top_loading_box_grouped_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-top-loading-box-grouped-{condition}.png',
		sig_top_loading_box_summary_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-top-loading-box-summary-{condition}.png',
		sig_top_count_box_grouped_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-top-count-box-grouped-{condition}.png',
		sig_top_count_box_summary_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-top-count-box-summary-{condition}.png',
		sig_top_freq_bar_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-top-freq-bar-{condition}.png',
	resources:
		mem_mb = 2000
	threads: 1
	script:
		'signature-analysis/visualize-liger-signature-loading-analysis.R'

rule liger_signature_loading_quantile_analysis: 
	input:
		sig_loading_mtx_long = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-long-form-{condition}.tsv',
	params:
	output:
		sig_loading_quantiles_cohort = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-cohort-{condition}.tsv',
		sig_loading_quantiles_sample = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-sample-{condition}.tsv',
		sig_loading_quantiles_signature = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-signature-{condition}.tsv',
		sig_loading_quantiles_signature_cohort = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-signature-cohort-{condition}.tsv',
		sig_loading_quantiles_signature_sample = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-signature-sample-{condition}.tsv',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'signature-analysis/liger-signature-loading-quantile-analysis.R'

rule visualize_liger_signature_loading_quantile_analysis: 
	input:
		sig_loading_quantiles_cohort = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-cohort-{condition}.tsv',
		sig_loading_quantiles_sample = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-sample-{condition}.tsv',
		sig_loading_quantiles_signature = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-signature-{condition}.tsv',
		sig_loading_quantiles_signature_cohort = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-signature-cohort-{condition}.tsv',
		sig_loading_quantiles_signature_sample = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-signature-sample-{condition}.tsv',
	params:
	output:
		sig_loading_quantiles_bar_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-single-var-{condition}.png',
		sig_loading_quantiles_box_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-multi-var-{condition}.png',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'signature-analysis/visualize-liger-signature-loading-quantile-analysis.R'

rule extract_liger_signature_loading_profiles:
	input:
		sig_top_freq = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-top-frequency-{condition}.tsv',
		sig_loading_patient_summary = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-patient-summary-{condition}.tsv',
		sig_activation_freq = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-activation-frequency-{condition}.tsv',
	params:
		
	output:
		sig_top_freq_profiles = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{condition}/{subtype}-signature-loading-profiles-top-frequency-{condition}.tsv',
		sig_loading_mean_profiles = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{condition}/{subtype}-signature-loading-profiles-loading-mean-{condition}.tsv',
		sig_loading_median_profiles = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{condition}/{subtype}-signature-loading-profiles-loading-median-{condition}.tsv',
		sig_act_freq_profiles = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{condition}/{subtype}-signature-loading-profiles-activation-frequency-{condition}.tsv',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'signature-analysis/extract-liger-signature-loading-profiles.R'

rule visualize_liger_signature_loading_profiles:
	input:
		sig_profiles = resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{condition}/{subtype}-signature-loading-profiles-{profile}-{condition}.tsv',
	params:
		
	output:
		sig_profiles_plot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{condition}/{subtype}-signature-loading-profiles-{profile}-{condition}.png',
		sig_profiles_corrplot = figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{condition}/profile-correlations/{subtype}-signature-loading-profiles-{profile}-correlation-{condition}.png',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'signature-analysis/visualize-liger-signature-loading-profiles.R'

