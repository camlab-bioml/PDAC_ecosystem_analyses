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
		patient_profiles_collapsed = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed.tsv',
		patient_profiles_collapsed_scored_val = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed-scored-validation/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-scored-validation.tsv',
	params:
	output:
		patient_profiles_collapsed_meta = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-metadata.tsv',
		patient_profiles_collapsed_cleaned_scaled = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-cleaned-scaled.tsv',
		patient_profiles_collapsed_scored_val_meta = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed-scored-validation/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-scored-validation-metadata.tsv',
		patient_profiles_collapsed_scored_val_cleaned_scaled = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed-scored-validation/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-scored-validation-cleaned-scaled.tsv',
		patient_profiles_correlation_data_frame = resultoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/patient-{compartment}-signature-profiles-{profile}-correlation-data-frame.tsv',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'patient-analysis/liger-signature-correlation-comparison.R'

rule visualize_liger_signature_correlation_comparison:
	input:
		celltype_pal = figureoutput + 'celltype-palette.rds',
		cohort_pal = figureoutput + 'cohort-palette.rds',
		cell_type_rename = config['figure1']['cell_type_rename_csv'],
		ambient_sigs = config['figure2']['ambient_sigs'],
		patient_profiles_collapsed_cleaned_scaled = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-cleaned-scaled.tsv',
		patient_profiles_collapsed_meta = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-metadata.tsv',
		patient_profiles_collapsed_scored_val_cleaned_scaled = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed-scored-validation/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-scored-validation-cleaned-scaled.tsv',
		patient_profiles_collapsed_scored_val_meta = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/collapsed-scored-validation/{profile}/patient-{compartment}-signature-profiles-{profile}-collapsed-scored-validation-metadata.tsv',
		patient_profiles_correlation_data_frame = resultoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/patient-{compartment}-signature-profiles-{profile}-correlation-data-frame.tsv',
	params:
	output:
		patient_profiles_collapsed_cleaned_sclaed_heatmap = figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/{compartment}-signature-profiles-{profile}-cleaned-scaled-heatmap.png',
		patient_profiles_collapsed_scored_val_cleaned_scaled_heatmap = figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/{compartment}-signature-profiles-{profile}-collapsed-scored-validation-cleaned-scaled-heatmap.png',
		patient_profiles_combined_cleaned_scaled_heatmap = figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/{compartment}-signature-profiles-{profile}-combined-cleaned-scaled-heatmap.png',
		patient_profiles_correlation_comparison_overall = figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/{compartment}-signature-profiles-{profile}-correlation-comparison-overall.png',
		patient_profiles_correlation_comparison_intra_celltype = figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/{compartment}-signature-profiles-{profile}-correlation-comparison-intra-celltype.png',
		patient_profiles_correlation_comparison_inter_celltype = figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/{compartment}-signature-profiles-{profile}-correlation-comparison-inter-celltype.png',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'patient-analysis/visualize-liger-signature-correlation-comparison.R'

rule stan_make_data:
	input:
		celltype_sig_loading_mtxs = lambda wildcards: expand(resultoutput + 'LIGER/signature-analysis/{celltype}/loading-matrices/{celltype}-signature-loading-' + wildcards.condition + '.tsv', celltype = list(stan_ct_config[wildcards.scope])),
		sig_loading_profiles = resultoutput + "LIGER/patient-analysis/patient-signature-profiles/{condition}/loading-mean/patient-full-signature-profiles-loading-mean-{condition}.tsv",
	params:
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.scope]),
		number_of_niches = lambda wildcards: wildcards.numniche,
		min_number_of_cells_per_sample = config['stan']['min_number_of_cells_per_sample'],
		max_number_of_cells_per_sample = config['stan']['max_number_of_cells_per_sample'],
		prep_fig_path = figureoutput + 'LIGER/patient-analysis/stan/prep-figures/{condition}/',
	output:
		stan_data = resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-data-{scope}-env-{condition}-num-niches-{numniche}.rds',
		df_sig_mean = resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-df-sig-mean-{scope}-env-{condition}-num-niches-{numniche}.rds',
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
		fit_optim = resultoutput + 'LIGER/patient-analysis/stan/results/model-fit/{condition}/stan-model-fit-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
	resources:
		mem_mb = 100000
	threads: 16
	script:
		'patient-analysis/stan-run-model.R'

rule stan_extract_parameters:
	input:
		fit_optim = resultoutput + 'LIGER/patient-analysis/stan/results/model-fit/{condition}/stan-model-fit-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
		stan_data = resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-data-{scope}-env-{condition}-num-niches-{numniche}.rds',
		df_sig_mean = resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-df-sig-mean-{scope}-env-{condition}-num-niches-{numniche}.rds',
	params:
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.scope]),
		number_of_niches = lambda wildcards: wildcards.numniche,
		nIter = lambda wildcards: wildcards.niter,
	output:
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
		patient_specific_modelled_mu = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/patient-specific-modelled-mu-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
		microenvironment_niche_factors = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/microenvironment-niche-factors-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
		niche_factor_loadings = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/niche-factor-loadings-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
		intrinsic_covariance_matrices = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/intrinsic-covariance-matrices-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds',
	params:
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.scope]),
		number_of_niches = lambda wildcards: wildcards.numniche,
		nIter = lambda wildcards: wildcards.niter,
	output:
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
		microenvironment_niche_factors_collapsed = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed/microenvironment-niche-factors-{scope}-env-collapsed-num-niches-{numniche}-niter-{niter}.rds',
		microenvironment_niche_factors_scored_val = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed-scored-validation/microenvironment-niche-factors-{scope}-env-collapsed-scored-validation-num-niches-{numniche}-niter-{niter}.rds',
		intrinsic_covariance_matrices_collapsed = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed/intrinsic-covariance-matrices-{scope}-env-collapsed-num-niches-{numniche}-niter-{niter}.rds',
		intrinsic_covariance_matrices_scored_val = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed-scored-validation/intrinsic-covariance-matrices-{scope}-env-collapsed-scored-validation-num-niches-{numniche}-niter-{niter}.rds',
	params:
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.scope]),
		number_of_niches = lambda wildcards: wildcards.numniche,
		nIter = lambda wildcards: wildcards.niter,
	output:
		microenvironment_niche_factors_plot = figureoutput + 'LIGER/patient-analysis/stan/results/parameter-comparison/microenvironment-niche-factors-{scope}-env-collapsed-and-scored-validation-num-niches-{numniche}-niter-{niter}-sig-loading.png',
		microenvironment_niche_factors_correlation_plot = figureoutput + 'LIGER/patient-analysis/stan/results/parameter-comparison/microenvironment-niche-factors-{scope}-env-collapsed-vs-scored-validation-num-niches-{numniche}-niter-{niter}-correlation.png',
		intrinsic_covariance_matrices_correlation_plot = figureoutput + 'LIGER/patient-analysis/stan/results/parameter-comparison/intrinsic-covariance-matrices-{scope}-env-collapsed-vs-scored-validation-num-niches-{numniche}-niter-{niter}-correlation.png',
		intrinsic_covariance_matrices_correlation_bar_plot = figureoutput + 'LIGER/patient-analysis/stan/results/parameter-comparison/intrinsic-covariance-matrices-{scope}-env-collapsed-vs-scored-validation-num-niches-{numniche}-niter-{niter}-correlation-bar.png',
	resources:
		mem_mb = 1000
	threads: 1
	script:
		'patient-analysis/stan-make-comparison-plots.R'