rule construct_liger_signature_patient_profiles:
	input:
		celltype_sig_profiles = lambda wildcards: expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/' + wildcards.condition + '/{subtype}-signature-loading-profiles-' + wildcards.profile + '-' + wildcards.condition + '.tsv', subtype = list(ct_for_patient_profile[wildcards.compartment])),
	params:
	output:
		patient_profiles = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}.tsv',
		patient_profiles_long = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-longform-{profile}-{condition}.tsv',
		patient_profiles_corr = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-correlation-{profile}-{condition}.tsv',
	resources:	
		mem_mb = 1000
	threads: 1
	script: 
		'patient-analysis/construct-liger-signature-patient-profiles.R'

rule visualize_liger_signature_patient_profiles:
	input:
		patient_profiles = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}.tsv',
		patient_profiles_long = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-longform-{profile}-{condition}.tsv',
		patient_profiles_corr = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-correlation-{profile}-{condition}.tsv',
	params:
		patient_profiles_individual_circlized_bar_plot_dir = figureoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}-individual-circlized-bar-plot/',
	output:
		patient_profiles_heatmap_grouped = figureoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}-heatmap-grouped.png',
		patient_profiles_heatmap_clustered = figureoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}-heatmap-clustered.png',
		patient_profiles_corrplot = figureoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}-corrplot.png',
		patient_profiles_stacked_bar_plot = figureoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}-stacked-bar-plot.png',
		patient_profiles_circlized_bar_plot = figureoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}-circlized-bar-plot.png',
	resources:	
		mem_mb = 1000
	threads: 1
	script: 
		'patient-analysis/visualize-liger-signature-patient-profiles.R'

rule liger_signature_collapsed_vs_collapsed_scored_validation:
	input:
		patient_profiles_corr_collapsed = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed/{profile}/patient-{compartment}-signature-profiles-correlation-{profile}-collapsed.tsv',
		patient_profiles_corr_collapsed_scored_val = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed-scored-validation/{profile}/patient-{compartment}-signature-profiles-correlation-{profile}-collapsed-scored-validation.tsv',
		patient_profiles_collapsed_scored_val = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed-scored-validation/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-scored-validation.tsv',
	params:
		compare_per_cohort = config['patient_profiles']['compare_collapsed_and_rescored_validation_per_cohort'],
		comparison_results_for_individual_cohorts_dir = resultoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-correlation-individual-validation-cohort-comparison/',
	output:
		sig_loading_corr_comparison = resultoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-correlation-comparison.tsv',
		sig_loading_corr_sign_comparison = resultoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-correlation-sign-comparison.tsv',
		sig_loading_corr_mean_comparison = resultoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-correlation-mean-comparison.tsv',
	resources: 
		mem_mb = 1000
	threads: 1
	script:
		'patient-analysis/liger-signature-collapsed-vs-collapsed-scored-validation.R'

rule visualize_liger_signature_collapsed_vs_collapsed_scored_validation:
	input:
		sig_loading_corr_comparison = resultoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-correlation-comparison.tsv',
		sig_loading_corr_sign_comparison = resultoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-correlation-sign-comparison.tsv',
		sig_loading_corr_mean_comparison = resultoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-correlation-mean-comparison.tsv',
	params:
		compare_per_cohort = config['patient_profiles']['compare_collapsed_and_rescored_validation_per_cohort'],
		comparison_results_for_individual_cohorts_dir = resultoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-correlation-individual-validation-cohort-comparison/',
		comparison_plots_for_individual_cohorts_dir = figureoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-correlation-individual-validation-cohort-comparison/',
	output:
		scatter_plot_all_intercell_sig_pairs = figureoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-all-intercell-sig-pair-correlation-comparison.png',
		scatter_plot_intercell_sig_pair_means = figureoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-intercell-sig-pair-correlation-mean-comparison.png',
		bar_plot_all_sig_pair_corr_sign_comparison = figureoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-all-sig-pair-correlation-sign-comparison.png',
		bar_plot_sig_pair_corr_same_sign_freq = figureoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-sig-pair-correlation-same-sign-frequency.png',
	resources: 
		mem_mb = 1000
	threads: 1
	script:
		'patient-analysis/visualize-liger-signature-collapsed-vs-collapsed-scored-validation.R'

rule liger_signature_correlation_comparison:
	input: 
		sig_interpretation = config['figure2']['sig_interpretation'],
		patient_profiles_collapsed = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed.tsv',
		patient_profiles_collapsed_scored_val = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed-scored-validation/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-scored-validation.tsv',
		validated_sig_df = lambda wildcards: expand(resultoutput + 'LIGER/signature-analysis/{celltype}/signature-validation/{celltype}-validated-signatures.tsv', celltype = list(ct_for_patient_profile[wildcards.compartment])),
		signature_collapse_guide = resultoutput + 'LIGER/signature-analysis/signature-collapse-guide.tsv',
	params:
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.compartment]),
	output:
		patient_profiles_collapsed_meta = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-metadata.tsv',
		patient_profiles_collapsed_cleaned_scaled = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-cleaned-scaled.tsv',
		patient_profiles_collapsed_scored_val_meta = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed-scored-validation/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-scored-validation-metadata.tsv',
		patient_profiles_collapsed_scored_val_cleaned_scaled = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed-scored-validation/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-scored-validation-cleaned-scaled.tsv',
		patient_profiles_correlation_data_frame = resultoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/patient-{compartment}-signature-profiles-{profile}-correlation-data-frame.tsv',
		signature_correlation_comparison_data_frame_list = resultoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/patient-{compartment}-signature-{profile}-correlation-comparison-data-frame-list.rds',
		signature_cooccurrence_agreement_data_frame = resultoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/patient-{compartment}-signature-{profile}-cooccurrence-agreement-data-frame.tsv',
	resources:
		mem_mb = 2000
	threads: 1
	script:
		'patient-analysis/liger-signature-correlation-comparison.R'

rule visualize_liger_signature_correlation_comparison:
	input:
		celltype_pal = figureoutput + 'celltype-palette.rds',
		cohort_pal = figureoutput + 'cohort-palette.rds',
		sig_interpretation = config['figure2']['sig_interpretation'],
		cell_type_rename = config['figure1']['cell_type_rename_csv'],
		ambient_sigs = config['figure2']['ambient_sigs'],
		patient_profiles_collapsed_cleaned_scaled = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-cleaned-scaled.tsv',
		patient_profiles_collapsed_meta = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-metadata.tsv',
		patient_profiles_collapsed_scored_val_cleaned_scaled = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed-scored-validation/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-scored-validation-cleaned-scaled.tsv',
		patient_profiles_collapsed_scored_val_meta = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed-scored-validation/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-scored-validation-metadata.tsv',
		patient_profiles_correlation_data_frame = resultoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/patient-{compartment}-signature-profiles-{profile}-correlation-data-frame.tsv',
		signature_correlation_comparison_data_frame_list = resultoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/patient-{compartment}-signature-{profile}-correlation-comparison-data-frame-list.rds',
		signature_cooccurrence_agreement_data_frame = resultoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/patient-{compartment}-signature-{profile}-cooccurrence-agreement-data-frame.tsv',
	params:
		sig_cooccurring_plot_dir = figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/sig-cooccur/',
	output:
		patient_profiles_collapsed_cleaned_sclaed_heatmap = figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/{compartment}-signature-profiles-{profile}-cleaned-scaled-heatmap.png',
		patient_profiles_collapsed_scored_val_cleaned_scaled_heatmap = figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/{compartment}-signature-profiles-{profile}-collapsed-scored-validation-cleaned-scaled-heatmap.png',
		patient_profiles_combined_cleaned_scaled_heatmap = figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/{compartment}-signature-profiles-{profile}-combined-cleaned-scaled-heatmap.png',
		patient_profiles_correlation_comparison_overall = figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/{compartment}-signature-profiles-{profile}-correlation-comparison-overall.png',
		patient_profiles_correlation_comparison_intra_celltype = figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/{compartment}-signature-profiles-{profile}-correlation-comparison-intra-celltype.png',
		patient_profiles_correlation_comparison_inter_celltype = figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/{compartment}-signature-profiles-{profile}-correlation-comparison-inter-celltype.png',
		signature_cooccurrence_agreement_plot = figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/{compartment}-signature-{profile}-cooccurrence-agreement-barplot.png',
	resources:
		mem_mb = 4000
	threads: 1
	script:
		'patient-analysis/visualize-liger-signature-correlation-comparison.R'

rule visualize_lige_signature_gene_expression_analysis:
	input:
		patient_profiles_correlation_data_frame = resultoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/patient-{compartment}-signature-profiles-{profile}-correlation-data-frame.tsv',
		sig_gene_expr_mtx_dis = resultoutput + 'LIGER/signature-analysis/{subtype}/gene-expression-analysis/{condition}/{subtype}-signature-gene-expression-collapsed.tsv',
		sig_gene_expr_mtx_val = resultoutput + 'LIGER/signature-analysis/{subtype}/gene-expression-analysis/{condition}/{subtype}-signature-gene-expression-collapsed-scored-validation.tsv',
	params:
	output:
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'patient-analysis/visualize-liger-signature-gene-expression-analysis.R'

rule figure_4_preparation:
	input: 
		patient_profiles_correlation_data_frame = resultoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/patient-{compartment}-signature-profiles-{profile}-correlation-data-frame.tsv',
		signature_correlation_comparison_data_frame_list = resultoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/patient-{compartment}-signature-{profile}-correlation-comparison-data-frame-list.rds',
		signature_cooccurrence_agreement_data_frame = resultoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/patient-{compartment}-signature-{profile}-cooccurrence-agreement-data-frame.tsv',
	params:
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.compartment]),
	output:
		patient_profiles_correlation_data_frame_full_and_intra = resultoutput + 'LIGER/patient-analysis/figure-4/{profile}/patient-{compartment}-signature-profiles-{profile}-correlation-data-frame-full-and-intra.tsv',
		patient_profiles_correlation_data_frame_inter = resultoutput + 'LIGER/patient-analysis/figure-4/{profile}/patient-{compartment}-signature-profiles-{profile}-correlation-data-frame-inter.tsv',
		signature_correlation_comparison_data_frame_list_examples = resultoutput + 'LIGER/patient-analysis/figure-4/{profile}/patient-{compartment}-signature-{profile}-correlation-comparison-data-frame-list-examples.rds',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'patient-analysis/figure-4-preparation.R'

rule draw_figure_4:
	input:
		celltype_pal = figureoutput + 'celltype-palette.rds',
		cohort_pal = figureoutput + 'cohort-palette.rds',
		sig_interpretation = config['figure2']['sig_interpretation'],
		cell_type_rename = config['figure1']['cell_type_rename_csv'],
		ambient_sigs = config['figure2']['ambient_sigs'],
		schematic = config['figure4']['schematic_plot'],
		patient_profiles_correlation_data_frame_full_and_intra = resultoutput + 'LIGER/patient-analysis/figure-4/loading-mean/patient-{compartment}-signature-profiles-loading-mean-correlation-data-frame-full-and-intra.tsv',
		#patient_profiles_correlation_data_frame_inter = resultoutput + 'LIGER/patient-analysis/figure-4/loading-mean/patient-{compartment}-signature-profiles-loading-mean-correlation-data-frame-inter.tsv',
		signature_correlation_comparison_data_frame_list_examples = resultoutput + 'LIGER/patient-analysis/figure-4/loading-mean/patient-{compartment}-signature-loading-mean-correlation-comparison-data-frame-list-examples.rds',
		signature_cooccurrence_agreement_data_frame = resultoutput + 'LIGER/patient-analysis/signature-correlation-comparison/loading-mean/patient-{compartment}-signature-loading-mean-cooccurrence-agreement-data-frame.tsv',
	params:
		figure4_a_width = config['figure4']['figure4_a_width'],
		figure4_a_height = config['figure4']['figure4_a_height'],
		figure4_b_width = config['figure4']['figure4_b_width'],
		figure4_b_height = config['figure4']['figure4_b_height'],
		figure4_c_width = config['figure4']['figure4_c_width'],
		figure4_c_height = config['figure4']['figure4_c_height'],
		figure4_d_width = config['figure4']['figure4_d_width'],
		figure4_d_height = config['figure4']['figure4_d_height'],
		figure4_e_width = config['figure4']['figure4_e_width'],
		figure4_e_height = config['figure4']['figure4_e_height'],
		figure4_f_width = config['figure4']['figure4_f_width'],
		figure4_f_height = config['figure4']['figure4_f_height'],
		figure4_width = config['figure4']['figure4_width'],
		figure4_height = config['figure4']['figure4_height'],
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.compartment]),
	output:
		figure4_a = figureoutput + 'LIGER/patient-analysis/figure-4/{compartment}/loading-mean/figure-4-a-loading-mean-{compartment}.png',
		figure4_b = figureoutput + 'LIGER/patient-analysis/figure-4/{compartment}/loading-mean/figure-4-b-loading-mean-{compartment}.png',
		figure4_c = figureoutput + 'LIGER/patient-analysis/figure-4/{compartment}/loading-mean/figure-4-c-loading-mean-{compartment}.png',
		figure4_d = figureoutput + 'LIGER/patient-analysis/figure-4/{compartment}/loading-mean/figure-4-d-loading-mean-{compartment}.png',
		figure4_e = figureoutput + 'LIGER/patient-analysis/figure-4/{compartment}/loading-mean/figure-4-e-loading-mean-{compartment}.png',
		figure4_f = figureoutput + 'LIGER/patient-analysis/figure-4/{compartment}/loading-mean/figure-4-f-loading-mean-{compartment}.png',
		figure4_png = figureoutput + 'LIGER/patient-analysis/figure-4/{compartment}/loading-mean/figure-4-loading-mean-{compartment}.png',
		figure4_pdf = figureoutput + 'LIGER/patient-analysis/figure-4/{compartment}/loading-mean/figure-4-loading-mean-{compartment}.pdf',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'patient-analysis/draw-figure-4.R'

rule draw_figure_2_new:
	input:
		celltype_pal = figureoutput + 'celltype-palette.rds',
		cohort_pal = figureoutput + 'cohort-palette.rds',
		sig_interpretation = config['figure2']['sig_interpretation'],
		cell_type_rename = config['figure1']['cell_type_rename_csv'],
		ambient_sigs = config['figure2']['ambient_sigs'],
		schematic = "resources/New_Cooccurrence_Schematic.pdf",
		patient_profiles_correlation_data_frame_full_and_intra = resultoutput + 'LIGER/patient-analysis/figure-4/loading-mean/patient-{compartment}-signature-profiles-loading-mean-correlation-data-frame-full-and-intra.tsv',
		#patient_profiles_correlation_data_frame_inter = resultoutput + 'LIGER/patient-analysis/figure-4/loading-mean/patient-{compartment}-signature-profiles-loading-mean-correlation-data-frame-inter.tsv',
		signature_correlation_comparison_data_frame_list_examples = resultoutput + 'LIGER/patient-analysis/figure-4/loading-mean/patient-{compartment}-signature-loading-mean-correlation-comparison-data-frame-list-examples.rds',
		signature_cooccurrence_agreement_data_frame = resultoutput + 'LIGER/patient-analysis/signature-correlation-comparison/loading-mean/patient-{compartment}-signature-loading-mean-cooccurrence-agreement-data-frame.tsv',
	params:
		figure2_a_width = config['figure4']['figure4_a_width'],
		figure2_a_height = config['figure4']['figure4_a_height'],
		figure2_b_width = 5.1,
		figure2_b_height = 5,
		figure2_c_width = config['figure4']['figure4_c_width'],
		figure2_c_height = config['figure4']['figure4_c_height'],
		figure2_d_width = config['figure4']['figure4_d_width'],
		figure2_d_height = config['figure4']['figure4_d_height'],
		figure2_e_width = config['figure4']['figure4_e_width'],
		figure2_e_height = config['figure4']['figure4_e_height'],
		figure2_f_width = config['figure4']['figure4_f_width'],
		figure2_f_height = config['figure4']['figure4_f_height'],
		figure2_width = 14,
		figure2_height = 16,
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.compartment]),
	output:
		figure2_a = figureoutput + 'LIGER/patient-analysis/figure-2-new/{compartment}/loading-mean/figure-2-a-loading-mean-{compartment}.png',
		figure2_b = figureoutput + 'LIGER/patient-analysis/figure-2-new/{compartment}/loading-mean/figure-2-b-loading-mean-{compartment}.png',
		figure2_c = figureoutput + 'LIGER/patient-analysis/figure-2-new/{compartment}/loading-mean/figure-2-c-loading-mean-{compartment}.png',
		figure2_d = figureoutput + 'LIGER/patient-analysis/figure-2-new/{compartment}/loading-mean/figure-2-d-loading-mean-{compartment}.png',
		figure2_e = figureoutput + 'LIGER/patient-analysis/figure-2-new/{compartment}/loading-mean/figure-2-e-loading-mean-{compartment}.png',
		figure2_f = figureoutput + 'LIGER/patient-analysis/figure-2-new/{compartment}/loading-mean/figure-2-f-loading-mean-{compartment}.png',
		figure2_png = figureoutput + 'LIGER/patient-analysis/figure-2-new/{compartment}/loading-mean/figure-2-loading-mean-{compartment}.png',
		figure2_pdf = figureoutput + 'LIGER/patient-analysis/figure-2-new/{compartment}/loading-mean/figure-2-loading-mean-{compartment}.pdf',
		signature_metrics_summary = resultoutput + 'LIGER/patient-analysis/figure-2-new/{compartment}/loading-mean/signature-metrics-summary.tsv',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'patient-analysis/draw-figure-2-new.R'


rule stan_make_data:
	input:
		celltype_sig_loading_mtxs = lambda wildcards: expand(resultoutput + 'LIGER/signature-analysis/{celltype}/loading-matrices/{celltype}-signature-loading-' + wildcards.condition + '.tsv', celltype = list(stan_ct_config[wildcards.scope])),
		sig_loading_profiles = resultoutput + "LIGER/patient-analysis/patient-signature-profiles/{condition}/loading-mean/patient-full-signature-profiles-loading-mean-{condition}.tsv",
		ambient_sigs = config['figure2']['ambient_sigs'],
	params:
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.scope]),
		number_of_niches = lambda wildcards: wildcards.numniche,
		min_number_of_cells_per_sample = config['stan']['min_number_of_cells_per_sample'],
		max_number_of_cells_per_sample = config['stan']['max_number_of_cells_per_sample'],
		prep_fig_path = figureoutput + 'LIGER/patient-analysis/stan/prep-figures/{condition}/',
	output:
		stan_data = resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-data-{scope}-env-{condition}-num-niches-{numniche}.rds',
		df_sig_mean = resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-df-sig-mean-{scope}-env-{condition}-num-niches-{numniche}.rds',
		samples_encodings = resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-samples-encodings-{scope}-env-{condition}-num-niches-{numniche}.tsv',
		sigs_encodings = resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-sigs-encodings-{scope}-env-{condition}-num-niches-{numniche}.tsv',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'patient-analysis/stan-make-data.R'

rule stan_write_model:
	input:
		stan_data = resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-data-{scope}-env-{condition}-num-niches-{numniche}.rds',
	params:
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.scope]),
		lambda_for_orthogonality = config['stan']['lambda_for_orthogonality'],
		stan_model_dir = resultoutput + 'LIGER/patient-analysis/stan/model/',
		stan_model_basename = 'model-{scope}-env-{condition}-num-niches-{numniche}-' + config['stan']['stan_model_basename'],
	output:
		stan_model = resultoutput + 'LIGER/patient-analysis/stan/model/model-{scope}-env-{condition}-num-niches-{numniche}-' + config['stan']['stan_model_basename'] + '.stan',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'patient-analysis/stan-write-model.R'

rule stan_run_model:
	input:
		stan_data = resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-data-{scope}-env-{condition}-num-niches-{numniche}.rds',
		df_sig_mean = resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-df-sig-mean-{scope}-env-{condition}-num-niches-{numniche}.rds',
		stan_model = resultoutput + 'LIGER/patient-analysis/stan/model/model-{scope}-env-{condition}-num-niches-{numniche}-' + config['stan']['stan_model_basename'] + '.stan',
	params:
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.scope]),
		number_of_niches = lambda wildcards: wildcards.numniche,
		nIter = lambda wildcards: wildcards.niter,
		nWarmup = config['stan']['nWarmup'],
		nChains = config['stan']['nChains'],
		treeDepth = config['stan']['treeDepth'],
	output:
		params_init_value_list = resultoutput + 'LIGER/patient-analysis/stan/results/model-fit/{condition}/stan-paparams-init-value-list-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
		fit_optim = resultoutput + 'LIGER/patient-analysis/stan/results/model-fit/{condition}/stan-model-fit-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
	resources:
		mem_mb = 100000
	threads: 32
	script:
		'patient-analysis/stan-run-model.R'

rule stan_extract_parameters:
	input:
		params_init_value_list = resultoutput + 'LIGER/patient-analysis/stan/results/model-fit/{condition}/stan-paparams-init-value-list-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
		fit_optim = resultoutput + 'LIGER/patient-analysis/stan/results/model-fit/{condition}/stan-model-fit-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
		stan_data = resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-data-{scope}-env-{condition}-num-niches-{numniche}.rds',
		df_sig_mean = resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-df-sig-mean-{scope}-env-{condition}-num-niches-{numniche}.rds',
		sigs_encodings = resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-sigs-encodings-{scope}-env-{condition}-num-niches-{numniche}.tsv',
		samples_encodings = resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-samples-encodings-{scope}-env-{condition}-num-niches-{numniche}.tsv',
	params:
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.scope]),
		number_of_niches = lambda wildcards: wildcards.numniche,
		nIter = lambda wildcards: wildcards.niter,
	output:
		microenvironment_niche_factors_init_value = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/microenvironment-niche-factors-init-value-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
		niche_factor_loadings_init_value = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/niche-factor-loadings-init-value-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
		patient_specific_modelled_mu = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/patient-specific-modelled-mu-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
		microenvironment_niche_factors = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/microenvironment-niche-factors-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
		niche_factor_loadings = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/niche-factor-loadings-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
		intrinsic_covariance_matrices = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/intrinsic-covariance-matrices-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
	resources:
		mem_mb = 2000
	threads: 1
	script:
		'patient-analysis/stan-extract-parameters.R'

rule stan_make_plots:
	input:
		celltype_pal = figureoutput + 'celltype-palette.rds',
		cohort_pal = figureoutput + 'cohort-palette.rds',
		cell_type_rename = config['figure1']['cell_type_rename_csv'],
		sig_interpretation = config['figure2']['sig_interpretation'],
		ambient_sigs = config['figure2']['ambient_sigs'],
		samples_encodings = resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-samples-encodings-{scope}-env-{condition}-num-niches-{numniche}.tsv',
		microenvironment_niche_factors_init_value = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/microenvironment-niche-factors-init-value-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
		niche_factor_loadings_init_value = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/niche-factor-loadings-init-value-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
		patient_specific_modelled_mu = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/patient-specific-modelled-mu-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
		microenvironment_niche_factors = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/microenvironment-niche-factors-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
		niche_factor_loadings = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/niche-factor-loadings-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
		intrinsic_covariance_matrices = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/intrinsic-covariance-matrices-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
	params:
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.scope]),
		number_of_niches = lambda wildcards: wildcards.numniche,
		nIter = lambda wildcards: wildcards.niter,
	output:
		microenvironment_niche_factors_init_value_plot = figureoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/microenvironment-niche-factors-init-value-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.png',
		niche_factor_loadings_init_value_plot = figureoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/niche-factor-loadings-init-value-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.png',
		patient_specific_modelled_mu_plot = figureoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/patient-specific-modelled-mu-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.png',
		microenvironment_niche_factors_plot = figureoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/microenvironment-niche-factors-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.png',
		niche_factor_loadings_plot = figureoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/niche-factor-loadings-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.png',
		intrinsic_covariance_matrices_plot = figureoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/intrinsic-covariance-matrices-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.png',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'patient-analysis/stan-make-plots.R'

rule stan_make_comparison_plots:
	input:
		celltype_pal = figureoutput + 'celltype-palette.rds',
		cohort_pal = figureoutput + 'cohort-palette.rds',
		cell_type_rename = config['figure1']['cell_type_rename_csv'],
		sig_interpretation = config['figure2']['sig_interpretation'],
		ambient_sigs = config['figure2']['ambient_sigs'],
		microenvironment_niche_factors_init_value_collapsed = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed/microenvironment-niche-factors-init-value-{scope}-env-collapsed-num-niches-{numniche}-niter-{niter}.rds',
		microenvironment_niche_factors_init_value_scored_val = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed-scored-validation/microenvironment-niche-factors-init-value-{scope}-env-collapsed-scored-validation-num-niches-{numniche}-niter-{niter}.rds',
		niche_factor_loadings_init_value_collapsed = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed/niche-factor-loadings-init-value-{scope}-env-collapsed-num-niches-{numniche}-niter-{niter}.rds',
		niche_factor_loadings_init_value_scored_val = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed-scored-validation/niche-factor-loadings-init-value-{scope}-env-collapsed-scored-validation-num-niches-{numniche}-niter-{niter}.rds',
		microenvironment_niche_factors_collapsed = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed/microenvironment-niche-factors-{scope}-env-collapsed-num-niches-{numniche}-niter-{niter}.rds',
		microenvironment_niche_factors_scored_val = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed-scored-validation/microenvironment-niche-factors-{scope}-env-collapsed-scored-validation-num-niches-{numniche}-niter-{niter}.rds',
		intrinsic_covariance_matrices_collapsed = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed/intrinsic-covariance-matrices-{scope}-env-collapsed-num-niches-{numniche}-niter-{niter}.rds',
		intrinsic_covariance_matrices_scored_val = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed-scored-validation/intrinsic-covariance-matrices-{scope}-env-collapsed-scored-validation-num-niches-{numniche}-niter-{niter}.rds',
	params:
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.scope]),
		number_of_niches = lambda wildcards: wildcards.numniche,
		nIter = lambda wildcards: wildcards.niter,
	output:
		microenvironment_niche_factors_init_value_plot = figureoutput + 'LIGER/patient-analysis/stan/results/parameter-comparison/microenvironment-niche-factors-init-value-{scope}-env-collapsed-and-scored-validation-num-niches-{numniche}-niter-{niter}-sig-loading.png',
		microenvironment_niche_factors_plot = figureoutput + 'LIGER/patient-analysis/stan/results/parameter-comparison/microenvironment-niche-factors-{scope}-env-collapsed-and-scored-validation-num-niches-{numniche}-niter-{niter}-sig-loading.png',
		microenvironment_niche_factors_init_value_correlation_plot = figureoutput + 'LIGER/patient-analysis/stan/results/parameter-comparison/microenvironment-niche-factors-init-value-{scope}-env-collapsed-vs-scored-validation-num-niches-{numniche}-niter-{niter}-correlation.png',
		microenvironment_niche_factors_correlation_plot = figureoutput + 'LIGER/patient-analysis/stan/results/parameter-comparison/microenvironment-niche-factors-{scope}-env-collapsed-vs-scored-validation-num-niches-{numniche}-niter-{niter}-correlation.png',
		intrinsic_covariance_matrices_correlation_plot = figureoutput + 'LIGER/patient-analysis/stan/results/parameter-comparison/intrinsic-covariance-matrices-{scope}-env-collapsed-vs-scored-validation-num-niches-{numniche}-niter-{niter}-correlation.png',
		intrinsic_covariance_matrices_correlation_bar_plot = figureoutput + 'LIGER/patient-analysis/stan/results/parameter-comparison/intrinsic-covariance-matrices-{scope}-env-collapsed-vs-scored-validation-num-niches-{numniche}-niter-{niter}-correlation-bar.png',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'patient-analysis/stan-make-comparison-plots.R'

rule prepare_Zhou_clinical_and_meta_data: 
	input:
		clinical_data = config['zhou']['clinical_data'],
		meta_data = config['zhou']['meta_data'],
	params:
		#cohort = config['zhou']['cohort'],
	output:
		clinical_data_cleaned = resultoutput + 'LIGER/patient-analysis/zhou/clinical-data-cleaned.tsv',
		meta_data_cleaned = resultoutput + 'LIGER/patient-analysis/zhou/meta-data-cleaned.tsv',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'patient-analysis/prepare-Zhou-clinical-and-meta-data.R'

rule figure_5_preparation:
	input:
		microenvironment_niche_factors_collapsed = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed/microenvironment-niche-factors-{scope}-env-collapsed-num-niches-4-niter-8000.rds',
		microenvironment_niche_factors_scored_val = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed-scored-validation/microenvironment-niche-factors-{scope}-env-collapsed-scored-validation-num-niches-4-niter-8000.rds',
		niche_factor_loadings_collapsed = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed/niche-factor-loadings-{scope}-env-collapsed-num-niches-4-niter-8000.rds',
		niche_factor_loadings_scored_val = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed-scored-validation/niche-factor-loadings-{scope}-env-collapsed-scored-validation-num-niches-4-niter-8000.rds',
		intrinsic_covariance_matrices_collapsed = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed/intrinsic-covariance-matrices-{scope}-env-collapsed-num-niches-4-niter-8000.rds',
		intrinsic_covariance_matrices_scored_val = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed-scored-validation/intrinsic-covariance-matrices-{scope}-env-collapsed-scored-validation-num-niches-4-niter-8000.rds',
	params:
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.scope]),
		cov_celltype_to_plot = config['figure5']['cov_celltype_to_plot'],
	output:
		microenvironment_niche_factors_combined = resultoutput + 'LIGER/patient-analysis/figure-5/microenvironment-niche-factors-{scope}-env-combined.tsv',
		#microenvironment_niche_factors_corr = resultoutput + 'LIGER/patient-analysis/figure-5/microenvironment-niche-factors-{scope}-env-corr.tsv',
		microenvironment_niche_factor_loadings_to_compare = resultoutput + 'LIGER/patient-analysis/figure-5/microenvironment-niche-factor-{scope}-env-loadings-to-compare.tsv',
		intrinsic_covariance_matrices_correlation_data = resultoutput + 'LIGER/patient-analysis/figure-5/intrinsic-covariance-matrices-{scope}-env-corr.tsv',
		intrinsic_covariance_matrices_list_to_plot = resultoutput + 'LIGER/patient-analysis/figure-5/intrinsic-covariance-matrices-{scope}-env.rds',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'patient-analysis/figure-5-preparation.R'

rule draw_figure_5:
	input:
		celltype_pal = figureoutput + 'celltype-palette.rds',
		cohort_pal = figureoutput + 'cohort-palette.rds',
		sig_interpretation = config['figure2']['sig_interpretation'],
		cell_type_rename = config['figure1']['cell_type_rename_csv'],
		ambient_sigs = config['figure2']['ambient_sigs'],
		schematic = config['figure5']['schematic_plot'],
		panel_d = config['figure5']['panel_d'],
		panel_e = config['figure5']['panel_e'],
		microenvironment_niche_factors_combined = resultoutput + 'LIGER/patient-analysis/figure-5/microenvironment-niche-factors-{scope}-env-combined.tsv',
		#microenvironment_niche_factors_corr = resultoutput + 'LIGER/patient-analysis/figure-5/microenvironment-niche-factors-{scope}-env-corr.tsv',
		microenvironment_niche_factor_loadings_to_compare = resultoutput + 'LIGER/patient-analysis/figure-5/microenvironment-niche-factor-{scope}-env-loadings-to-compare.tsv',
		intrinsic_covariance_matrices_correlation_data = resultoutput + 'LIGER/patient-analysis/figure-5/intrinsic-covariance-matrices-{scope}-env-corr.tsv',
		intrinsic_covariance_matrices_list_to_plot = resultoutput + 'LIGER/patient-analysis/figure-5/intrinsic-covariance-matrices-{scope}-env.rds',
	params:
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.scope]),
		figure5_a_width = config['figure5']['figure5_a_width'],
		figure5_a_height = config['figure5']['figure5_a_height'],
		figure5_b_width = config['figure5']['figure5_b_width'],
		figure5_b_height = config['figure5']['figure5_b_height'],
		figure5_c_width = config['figure5']['figure5_c_width'],
		figure5_c_height = config['figure5']['figure5_c_height'],
		figure5_d_width = config['figure5']['figure5_d_width'],
		figure5_d_height = config['figure5']['figure5_d_height'],
		figure5_e_width = config['figure5']['figure5_e_width'],
		figure5_e_height = config['figure5']['figure5_e_height'],
		figure5_width = config['figure5']['figure5_width'],
		figure5_height = config['figure5']['figure5_height'],
	output:
		figure5_a = figureoutput + 'LIGER/patient-analysis/figure-5/{scope}/figure-5-a-{scope}.png',
		figure5_b = figureoutput + 'LIGER/patient-analysis/figure-5/{scope}/figure-5-b-{scope}.png',
		figure5_c = figureoutput + 'LIGER/patient-analysis/figure-5/{scope}/figure-5-c-{scope}.png',
		figure5_d = figureoutput + 'LIGER/patient-analysis/figure-5/{scope}/figure-5-d-{scope}.png',
		figure5_e = figureoutput + 'LIGER/patient-analysis/figure-5/{scope}/figure-5-e-{scope}.png',
		figure5_png = figureoutput + 'LIGER/patient-analysis/figure-5/{scope}/figure-5-{scope}.png',
		figure5_pdf = figureoutput + 'LIGER/patient-analysis/figure-5/{scope}/figure-5-{scope}.pdf',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'patient-analysis/draw-figure-5.R'

rule draw_figure_3_new:
	input:
		celltype_pal = figureoutput + 'celltype-palette.rds',
		cohort_pal = figureoutput + 'cohort-palette.rds',
		sig_interpretation = config['figure2']['sig_interpretation'],
		cell_type_rename = config['figure1']['cell_type_rename_csv'],
		ambient_sigs = config['figure2']['ambient_sigs'],
		schematic = 'resources/Visium_validation_schematic.pdf',
		panel_b = 'resources/ecosystem_disval_corr.pdf',
		panel_c = 'resources/moranI_eyeplot.pdf',
		panel_d = 'resources/cNMF_eyeplot.pdf',
		panel_e1 = 'resources/niches_1-HT242P1H1-S1Fc1U1.pdf',
		panel_e2 = 'resources/niches_B1_HT270P1-S1H1Fs5U1.pdf',
		panel_e3 = 'resources/niches_B1-HT264P1-S1H2Fc2U1.pdf',
		panel_e4 = 'resources/niches_D1-HT306P1-S1H1Fc2U1.pdf',
		panel_e = 'resources/Visium_spots_samples.pdf',
		panel_d_supp = config['figure5']['panel_d'],
		panel_e_supp = config['figure5']['panel_e'],
		microenvironment_niche_factors_combined = resultoutput + 'LIGER/patient-analysis/figure-5/microenvironment-niche-factors-{scope}-env-combined.tsv',
		#microenvironment_niche_factors_corr = resultoutput + 'LIGER/patient-analysis/figure-5/microenvironment-niche-factors-{scope}-env-corr.tsv',
		microenvironment_niche_factor_loadings_to_compare = resultoutput + 'LIGER/patient-analysis/figure-5/microenvironment-niche-factor-{scope}-env-loadings-to-compare.tsv',
		intrinsic_covariance_matrices_correlation_data = resultoutput + 'LIGER/patient-analysis/figure-5/intrinsic-covariance-matrices-{scope}-env-corr.tsv',
		intrinsic_covariance_matrices_list_to_plot = resultoutput + 'LIGER/patient-analysis/figure-5/intrinsic-covariance-matrices-{scope}-env.rds',
	params:
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.scope]),
		figure3_a_width = 5,
		figure3_a_height = 5,
		figure3_b_width = config['figure5']['figure5_b_width'],
		figure3_b_height = config['figure5']['figure5_b_height'],
		figure3_c_width = config['figure5']['figure5_c_width'],
		figure3_c_height = config['figure5']['figure5_c_height'],
		figure3_d_width = config['figure5']['figure5_d_width'],
		figure3_d_height = config['figure5']['figure5_d_height'],
		figure3_e_width = config['figure5']['figure5_e_width'],
		figure3_e_height = config['figure5']['figure5_e_height'],
		figure3_width = 12,
		figure3_height = 14,
	output:
		figure3_a = figureoutput + 'LIGER/patient-analysis/figure-3-new/{scope}/figure-3-a-{scope}.png',
		figure3_b = figureoutput + 'LIGER/patient-analysis/figure-3-new/{scope}/figure-3-b-{scope}.png',
		figure3_c = figureoutput + 'LIGER/patient-analysis/figure-3-new/{scope}/figure-3-c-{scope}.png',
		figure3_d = figureoutput + 'LIGER/patient-analysis/figure-3-new/{scope}/figure-3-d-{scope}.png',
		figure3_e = figureoutput + 'LIGER/patient-analysis/figure-3-new/{scope}/figure-3-e-{scope}.png',
		figure3_c_supp = figureoutput + 'LIGER/patient-analysis/figure-3-new/{scope}/figure-3-c-supp-{scope}.pdf',
		figure3_d_supp = figureoutput + 'LIGER/patient-analysis/figure-3-new/{scope}/figure-3-d-supp-{scope}.png',
		figure3_e_supp = figureoutput + 'LIGER/patient-analysis/figure-3-new/{scope}/figure-3-e-supp-{scope}.png',
		figure3_png = figureoutput + 'LIGER/patient-analysis/figure-3-new/{scope}/figure-3-{scope}.png',
		figure3_pdf = figureoutput + 'LIGER/patient-analysis/figure-3-new/{scope}/figure-3-{scope}.pdf',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'patient-analysis/draw-figure-3-new.R'