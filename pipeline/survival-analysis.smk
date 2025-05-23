rule tidy_up_clinical_data:
	input:
		tcga_clin_data = config["survival"]["tcga_clin_data"],
		tcga_purist_subtype_calls = config["survival"]["tcga_purist_subtype_calls"],
		compass_clin_data = config["survival"]["compass_clin_data"],
		toronto_bulk_data = config["survival"]["toronto_bulk_data"],
		toronto_subtype_calls = config["survival"]["toronto_subtype_calls"],
		toronto_meta_data = config["survival"]["toronto_meta_data"],
		pancurx_ihc_data = config["survival"]["pancurx_ihc_data"],
	params:
	output:
		tidy_paad_clin_data = resultoutput + 'LIGER/survival-analysis/clinical-data/paad-clin-data.tsv',
		tidy_compass_clin_data = resultoutput + 'LIGER/survival-analysis/clinical-data/compass-clin-data.tsv',
		tidy_pancurx_clin_data = resultoutput + 'LIGER/survival-analysis/clinical-data/pancurx-clin-data.tsv',
		toronto_compass_clin_data_comparison = resultoutput + 'LIGER/survival-analysis/clinical-data/toronto-compass-clin-data-comparison.tsv',
	resources:	
		mem_mb = 10000
	threads: 1
	script: 
		'survival-analysis/tidy-up-clinical-data.R'

rule tidy_up_bulk_rna_data:
	input:
		paad_mrna_data = config["survival"]["paad_mrna_data"],
		paad_sample_list = config["survival"]["paad_sample_list"],
		toronto_bulk_data = config["survival"]["toronto_bulk_data"],
		tidy_compass_clin_data = resultoutput + 'LIGER/survival-analysis/clinical-data/compass-clin-data.tsv',
		tidy_pancurx_clin_data = resultoutput + 'LIGER/survival-analysis/clinical-data/pancurx-clin-data.tsv',
	params:
	output:
		tidy_paad_mrna_data = resultoutput + 'LIGER/survival-analysis/bulk-rna-data/paad-mrna-data.rds',
		tidy_pancurx_mrna_data = resultoutput + 'LIGER/survival-analysis/bulk-rna-data/pancurx-mrna-data.rds',
		tidy_compass_mrna_data = resultoutput + 'LIGER/survival-analysis/bulk-rna-data/compass-mrna-data.rds',
		tidy_compass_pa_mrna_data = resultoutput + 'LIGER/survival-analysis/bulk-rna-data/compass-pa-mrna-data.rds',
		tidy_compass_lv_mrna_data = resultoutput + 'LIGER/survival-analysis/bulk-rna-data/compass-lv-mrna-data.rds',
	resources:
		mem_mb = 10000
	threads: 1
	script: 
		'survival-analysis/tidy-up-bulk-rna-data.R'

rule signature_survival_analysis:
	input:
		cell_type_rename = config['figure1']['cell_type_rename_csv'],
		sig_interpretation = config['figure2']['sig_interpretation'],
		tidy_paad_mrna_data = resultoutput + 'LIGER/survival-analysis/bulk-rna-data/paad-mrna-data.rds',
		tidy_pancurx_mrna_data = resultoutput + 'LIGER/survival-analysis/bulk-rna-data/pancurx-mrna-data.rds',
		tidy_paad_clin_data = resultoutput + 'LIGER/survival-analysis/clinical-data/paad-clin-data.tsv',
		tidy_pancurx_clin_data = resultoutput + 'LIGER/survival-analysis/clinical-data/pancurx-clin-data.tsv',
		celltype_gene_loading_mtxs = lambda wildcards: expand(resultoutput + 'LIGER/signature-analysis/{celltype}/loading-matrices/{celltype}-gene-loading-collapsed.tsv', celltype = list(stan_ct_config[wildcards.scope])),
	params:
		celltypes = lambda wildcards: list(stan_ct_config[wildcards.scope]),
		score_used_for_clinical_association = config['survival']['score_used_for_clinical_association'],
	output:
		gsva_es_paad = resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/gsva-es-paad.rds',
		clin_data_for_plotting_paad = resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/clin-data-for-plotting-paad.rds',
		cm_curve_plot_list_paad = resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/cm-curve-plot-list-paad.rds',
		coxph_summary_paad = resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/coxph-summary-paad.rds',
		clin_assoc_paad = resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/clin-assoc-paad.rds',
		clin_assoc_paad_pval = resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/clin-assoc-paad-pval.rds',
		gsva_es_pancurx = resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/gsva-es-pancurx.rds',
		clin_data_for_plotting_pancurx = resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/clin-data-for-plotting-pancurx.rds',
		cm_curve_plot_list_pancurx = resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/cm-curve-plot-list-pancurx.rds',
		coxph_summary_pancurx = resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/coxph-summary-pancurx.rds',
		clin_assoc_pancurx = resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/clin-assoc-pancurx.rds',
		clin_assoc_pancurx_pval = resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/clin-assoc-pancurx-pval.rds',
	resources:
		mem_mb = 80000
	threads: 1
	conda: 
		"/home/campbell/cyu/miniconda3/envs/pdactme_survival",
	script: 
		'survival-analysis/signature-survival-analysis.R'

rule niche_survival_analysis:
	input:
		niche_factor_loadings = resultoutput + 'LIGER/patient-analysis/stan/results/parameters/collapsed/microenvironment-niche-factors-{scope}-env-collapsed-num-niches-4-niter-8000.rds',
		gsva_es_paad = resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/gsva-es-paad.rds',
		gsva_es_pancurx = resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/gsva-es-pancurx.rds',
		clin_data_for_plotting_paad = resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/clin-data-for-plotting-paad.rds',
		clin_data_for_plotting_pancurx = resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/clin-data-for-plotting-pancurx.rds',
	params:
		score_used_for_clinical_association = config['survival']['score_used_for_clinical_association'],
	output:
		gsva_es_paad_niche = resultoutput + 'LIGER/survival-analysis/{scope}/niche-survival-analysis/gsva-es-paad-niche.rds',
		sig_exprs_mean_paad_niche = resultoutput + 'LIGER/survival-analysis/{scope}/niche-survival-analysis/sig-exprs-mean-paad-niche.rds',
		clin_data_for_plotting_paad_niche = resultoutput + 'LIGER/survival-analysis/{scope}/niche-survival-analysis/clin-data-for-plotting-paad-niche.rds',
		cm_curve_plot_list_paad_niche = resultoutput + 'LIGER/survival-analysis/{scope}/niche-survival-analysis/cm-curve-plot-list-paad-niche.rds',
		coxph_summary_paad_niche = resultoutput + 'LIGER/survival-analysis/{scope}/niche-survival-analysis/coxph-summary-paad-niche.rds',
		clin_assoc_paad_niche = resultoutput + 'LIGER/survival-analysis/{scope}/niche-survival-analysis/clin-assoc-paad-niche.rds',
		clin_assoc_paad_pval_niche = resultoutput + 'LIGER/survival-analysis/{scope}/niche-survival-analysis/clin-assoc-paad-pval-niche.rds',
		gsva_es_pancurx_niche = resultoutput + 'LIGER/survival-analysis/{scope}/niche-survival-analysis/gsva-es-pancurx-niche.rds',
		sig_exprs_mean_pancurx_niche = resultoutput + 'LIGER/survival-analysis/{scope}/niche-survival-analysis/sig-exprs-mean-pancurx-niche.rds',
		clin_data_for_plotting_pancurx_niche = resultoutput + 'LIGER/survival-analysis/{scope}/niche-survival-analysis/clin-data-for-plotting-pancurx-niche.rds',
		cm_curve_plot_list_pancurx_niche = resultoutput + 'LIGER/survival-analysis/{scope}/niche-survival-analysis/cm-curve-plot-list-pancurx-niche.rds',
		coxph_summary_pancurx_niche = resultoutput + 'LIGER/survival-analysis/{scope}/niche-survival-analysis/coxph-summary-pancurx-niche.rds',
		clin_assoc_pancurx_niche = resultoutput + 'LIGER/survival-analysis/{scope}/niche-survival-analysis/clin-assoc-pancurx-niche.rds',
		clin_assoc_pancurx_pval_niche = resultoutput + 'LIGER/survival-analysis/{scope}/niche-survival-analysis/clin-assoc-pancurx-pval-niche.rds',
	resources:
		mem_mb = 1000
	threads: 1
	conda: 
		"/home/campbell/cyu/miniconda3/envs/pdactme_survival",
	script: 
		'survival-analysis/niche-survival-analysis.R'

rule visualize_niche_loadings_in_bulk_sample:
	input:
		gsva_es_paad_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/gsva-es-paad-niche.rds',
		sig_exprs_mean_paad_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/sig-exprs-mean-paad-niche.rds',
		gsva_es_pancurx_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/gsva-es-pancurx-niche.rds',
		sig_exprs_mean_pancurx_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/sig-exprs-mean-pancurx-niche.rds',
		tidy_pancurx_clin_data = resultoutput + 'LIGER/survival-analysis/clinical-data/pancurx-clin-data.tsv',
		tidy_paad_clin_data = resultoutput + 'LIGER/survival-analysis/clinical-data/paad-clin-data.tsv',
	params:
	output:
		figure_niche_loadings_paad = figureoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/niche-loadings-paad.png',
		figure_niche_loadings_pancurx = figureoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/niche-loadings-pancurx.png',
		figure_scatter_niche_loadings_paad = figureoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/scatter-niche-loadings-paad.png',
		figure_scatter_niche_loadings_pancurx = figureoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/scatter-niche-loadings-pancurx.png',
	resources:
		mem_mb = 1000
	threads: 1
	conda: 
		"/home/campbell/cyu/miniconda3/envs/pdactme_survival",
	script:
		'survival-analysis/visualize-niche-loadings-in-bulk-sample.R'

rule draw_figure_6:
	input:
		celltype_pal = figureoutput + 'celltype-palette.rds',
		cell_type_rename = config['figure1']['cell_type_rename_csv'],
		sig_interpretation = config['figure2']['sig_interpretation'],
		ambient_sigs = config['figure2']['ambient_sigs'],
		clin_data_for_plotting_paad = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/clin-data-for-plotting-paad.rds',
		clin_data_for_plotting_pancurx = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/clin-data-for-plotting-pancurx.rds',
		clin_data_for_plotting_paad_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/clin-data-for-plotting-paad-niche.rds',
		clin_data_for_plotting_pancurx_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/clin-data-for-plotting-pancurx-niche.rds',
		cm_curve_plot_list_paad = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/cm-curve-plot-list-paad.rds',
		cm_curve_plot_list_pancurx = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/cm-curve-plot-list-pancurx.rds',
		cm_curve_plot_list_paad_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/cm-curve-plot-list-paad-niche.rds',
		cm_curve_plot_list_pancurx_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/cm-curve-plot-list-pancurx-niche.rds',
		coxph_summary_paad = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/coxph-summary-paad.rds',
		coxph_summary_pancurx = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/coxph-summary-pancurx.rds',
		coxph_summary_paad_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/coxph-summary-paad-niche.rds',
		coxph_summary_pancurx_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/coxph-summary-pancurx-niche.rds',
		clin_assoc_paad = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/clin-assoc-paad.rds',
		clin_assoc_pancurx = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/clin-assoc-pancurx.rds',
		clin_assoc_paad_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/clin-assoc-paad-niche.rds',
		clin_assoc_pancurx_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/clin-assoc-pancurx-niche.rds',
		clin_assoc_paad_pval = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/clin-assoc-paad-pval.rds',
		clin_assoc_pancurx_pval = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/clin-assoc-pancurx-pval.rds',
		clin_assoc_paad_pval_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/clin-assoc-paad-pval-niche.rds',
		clin_assoc_pancurx_pval_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/clin-assoc-pancurx-pval-niche.rds',
	params:
		figure6_a_width = config['figure6']['figure6_a_width'],
		figure6_a_height = config['figure6']['figure6_a_height'],
		figure6_b_width = config['figure6']['figure6_b_width'],
		figure6_b_height = config['figure6']['figure6_b_height'],
		figure6_c_width = config['figure6']['figure6_c_width'],
		figure6_c_height = config['figure6']['figure6_c_height'],
		figure6_d_width = 11,
		figure6_d_height = 5,
		figure6_e_width = 15,
		figure6_e_height = 20,
		figure6_width = config['figure6']['figure6_width'],
		figure6_height = config['figure6']['figure6_height'],
	output:
		figure6_a = figureoutput + 'LIGER/survival-analysis/figure-6/figure-6-a.png',
		figure6_b = figureoutput + 'LIGER/survival-analysis/figure-6/figure-6-b.png',
		figure6_c = figureoutput + 'LIGER/survival-analysis/figure-6/figure-6-c.png',
		figure6_d = figureoutput + 'LIGER/survival-analysis/figure-6/figure-6-d.png',
		figure6_e = figureoutput + 'LIGER/survival-analysis/figure-6/figure-6-e.png',
		figure6_png = figureoutput + 'LIGER/survival-analysis/figure-6/figure-6.png',
		figure6_pdf = figureoutput + 'LIGER/survival-analysis/figure-6/figure-6.pdf',
	resources:
		mem_mb = 10000
	threads: 1
	conda: 
		"/home/campbell/cyu/miniconda3/envs/pdactme_survival",
	script: 
		'survival-analysis/draw-figure-6.R'

rule draw_figure_4_new:
	input:
		celltype_pal = figureoutput + 'celltype-palette.rds',
		cell_type_rename = config['figure1']['cell_type_rename_csv'],
		sig_interpretation = config['figure2']['sig_interpretation'],
		ambient_sigs = config['figure2']['ambient_sigs'],
		panel_a = config['figure5']['panel_e'],
		panel_b = figureoutput + 'LIGER/patient-analysis/figure-3-new/full/figure-3-c-supp-full.pdf',
		clin_data_for_plotting_paad = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/clin-data-for-plotting-paad.rds',
		clin_data_for_plotting_pancurx = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/clin-data-for-plotting-pancurx.rds',
		clin_data_for_plotting_paad_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/clin-data-for-plotting-paad-niche.rds',
		clin_data_for_plotting_pancurx_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/clin-data-for-plotting-pancurx-niche.rds',
		cm_curve_plot_list_paad = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/cm-curve-plot-list-paad.rds',
		cm_curve_plot_list_pancurx = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/cm-curve-plot-list-pancurx.rds',
		cm_curve_plot_list_paad_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/cm-curve-plot-list-paad-niche.rds',
		cm_curve_plot_list_pancurx_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/cm-curve-plot-list-pancurx-niche.rds',
		coxph_summary_paad = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/coxph-summary-paad.rds',
		coxph_summary_pancurx = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/coxph-summary-pancurx.rds',
		coxph_summary_paad_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/coxph-summary-paad-niche.rds',
		coxph_summary_pancurx_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/coxph-summary-pancurx-niche.rds',
		clin_assoc_paad = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/clin-assoc-paad.rds',
		clin_assoc_pancurx = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/clin-assoc-pancurx.rds',
		clin_assoc_paad_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/clin-assoc-paad-niche.rds',
		clin_assoc_pancurx_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/clin-assoc-pancurx-niche.rds',
		clin_assoc_paad_pval = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/clin-assoc-paad-pval.rds',
		clin_assoc_pancurx_pval = resultoutput + 'LIGER/survival-analysis/full/signature-survival-analysis/clin-assoc-pancurx-pval.rds',
		clin_assoc_paad_pval_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/clin-assoc-paad-pval-niche.rds',
		clin_assoc_pancurx_pval_niche = resultoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/clin-assoc-pancurx-pval-niche.rds',
	params:
		figure4_a_width = 8,
		figure4_a_height = 5,
		figure4_b_width = 11,
		figure4_b_height = 5,
		figure4_c_width = 6,
		figure4_c_height = 5,
		figure4_d_width = 11,
		figure4_d_height = 5,
		figure4_e_width = 10,
		figure4_e_height = 8,
		figure4_width = 11,
		figure4_height = 15,
	output:
		figure4_a = figureoutput + 'LIGER/survival-analysis/figure-4-new/figure-4-a.png',
		figure4_a_supp = figureoutput + 'LIGER/survival-analysis/figure-4-new/figure-4-a-supp.png',
		figure4_b = figureoutput + 'LIGER/survival-analysis/figure-4-new/figure-4-b.png',
		figure4_b_supp = figureoutput + 'LIGER/survival-analysis/figure-4-new/figure-4-b-supp.png',
		figure4_c = figureoutput + 'LIGER/survival-analysis/figure-4-new/figure-4-c.png',
		figure4_d = figureoutput + 'LIGER/survival-analysis/figure-4-new/figure-4-d.png',
		figure4_e_chan = figureoutput + 'LIGER/survival-analysis/figure-4-new/figure-4-e-chan.png',
		figure4_e_chan_data = resultoutput + 'LIGER/survival-analysis/figure-4-new/figure-4-e-chan-data.tsv',
		figure4_e_tcga = figureoutput + 'LIGER/survival-analysis/figure-4-new/figure-4-e-tcga.png',
		figure4_e_tcga_data = resultoutput + 'LIGER/survival-analysis/figure-4-new/figure-4-e-tcga-data.tsv',
		figure4_e_supp = figureoutput + 'LIGER/survival-analysis/figure-4-new/figure-4-e-supp.png',
		figure4_bar_supp_tcga = figureoutput + 'LIGER/survival-analysis/figure-4-new/figure-4-bar-supp_tcga.png',
		figure4_bar_supp_pancurx = figureoutput + 'LIGER/survival-analysis/figure-4-new/figure-4-bar-supp_pancurx.png',
		figure4_png = figureoutput + 'LIGER/survival-analysis/figure-4-new/figure-4.png',
		figure4_pdf = figureoutput + 'LIGER/survival-analysis/figure-4-new/figure-4.pdf',
	resources:
		mem_mb = 10000
	threads: 1
	conda: 
		"/home/campbell/cyu/miniconda3/envs/pdactme_survival",
	script: 
		'survival-analysis/draw-figure-4-new.R'