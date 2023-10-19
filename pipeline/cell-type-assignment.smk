rule singler_annotation_all_cells:
	input:
		sce_query = dataoutput + 'process-sce/scRNASeq-filtered-sce-{cohort}.rds',
		sce_ref = config['singler']['pdac_reference'],
	params:
		ref_label_field = config['singler']['pdac_ref_label_field'],
	output:
		sce = dataoutput + 'SingleR-annotation-all-cells/scRNASeq-SingleR-annotated-sce-{cohort}.rds',
		dframe = resultoutput + 'SingleR-annotation-all-cells/DFrame/cell-type-prediction-{cohort}.rds',
	resources:
		mem_mb = 20000
	threads: 12
	script:
		'cell-type-assignment/singler-annotation-all-cells.R'

rule visualize_singler_output_all_cells:
	input:
		sce = dataoutput + 'SingleR-annotation-all-cells/scRNASeq-SingleR-annotated-sce-{cohort}.rds',
		dframe = resultoutput + 'SingleR-annotation-all-cells/DFrame/cell-type-prediction-{cohort}.rds',
	params:
		ncol_violin = config['singler']['ncol_for_violin_plots'],
		celltype_label_field = config['singler']['celltype_label_field'],
		individual_celltype_marker_expression_heatmap_dir = figureoutput + 'SingleR-annotation-all-cells/celltype-marker-expression/{cohort}/',
	output:
		heatmap_assignment_score = figureoutput + 'SingleR-annotation-all-cells/celltype-assignment-score/heatmap-celltype-assignment-score-{cohort}.png',
		violin_delta_distribution = figureoutput + 'SingleR-annotation-all-cells/score-distribution/violin-celltype-delta-distribution-{cohort}.png',
		violin_score_distribution = figureoutput + 'SingleR-annotation-all-cells/score-distribution/violin-celltype-score-distribution-{cohort}.png',
	resources:
		mem_mb = 6000
	threads: 2
	script:
		'cell-type-assignment/visualize-singler-output-all-cells.R'

rule visualize_singler_output_dimred_all_cells:
	input:
		sce = dataoutput + 'SingleR-annotation-all-cells/scRNASeq-SingleR-annotated-sce-{cohort}.rds',
		dframe = resultoutput + 'SingleR-annotation-all-cells/DFrame/cell-type-prediction-{cohort}.rds',
	params:
		celltype_label_field = config['singler']['celltype_label_field'],
		dimred_to_plot = '{dimred}',
		metadata_to_color = '{metadata}',
		sce_assay_to_plot = config['singler']['sce_assay_to_plot'],
		individual_celltype_score_dimred_dir = figureoutput + 'SingleR-annotation-all-cells/celltype-assignment-score/individual-dimred/{cohort}/',
	output:
		dimred_cell_type = figureoutput+ 'SingleR-annotation-all-cells/annotated-dimred/{dimred}-cell-type-and-{metadata}-{cohort}.png',
	resources: 
		mem_mb = 6000
	threads: 2
	script:
		'cell-type-assignment/visualize-singler-output-dimred-all-cells.R'

rule plot_cell_type_markers_all_cells:
	input: 
		sce = dataoutput + 'SingleR-annotation-all-cells/scRNASeq-SingleR-annotated-sce-{cohort}.rds',
		dframe = resultoutput + 'SingleR-annotation-all-cells/DFrame/cell-type-prediction-{cohort}.rds',
	params: 
		dimred_to_plot = config['singler']['redim_to_plot'],
		celltype_label_field = config['singler']['celltype_label_field'],
		sce_assay_to_plot = config['singler']['sce_assay_to_plot'],
		marker = '{marker}',
		individual_plots_dir = figureoutput + 'SingleR-annotation-all-cells/celltype-markers/{cohort}/individual-plots/{marker}/', 
	output: 
		expression_plot = figureoutput + 'SingleR-annotation-all-cells/celltype-markers/{cohort}/{marker}.png',
	resources: 
		mem_mb = 3000
	threads: 1
	script: 
		'cell-type-assignment/plot-cell-type-markers-all-cells.R'

rule infercnv_preparation:
	input:
		sce_tumor = lambda wildcards: get_sce_dir(wildcards.celltype, dataoutput, celltypes, ct_zoomout, celltype_hierarchy) + 'scRNASeq-' + wildcards.celltype + '-sce-{cohort}.rds',
		sce_normal = lambda wildcards: get_sce_dir(wildcards.celltype, dataoutput, celltypes, ct_zoomout, celltype_hierarchy) + 'scRNASeq-' + wildcards.celltype + '-sce-{normcohort}.rds',
	params:
		celltype_label_field = config['singler']['celltype_label_field'],
	output:
		sce_combined = dataoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/scRNASeq-infercnv-prepared-sce-{celltype}-{cohort}-ref-{normcohort}.rds',
		sample_list = resultoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/sample-list.csv',
	resources: 
		mem_mb = 25000
	threads: 2
	script:
		'cell-type-assignment/infercnv-preparation.R'

rule run_infercnv:
	input:
		sce_combined = dataoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/scRNASeq-infercnv-prepared-sce-{celltype}-{cohort}-ref-{normcohort}.rds',
		gene_order_file = config['infercnv']['gene_order_file'],
	params:
		celltype_label_field = config['singler']['celltype_label_field'],
		infercnv_outdir = resultoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/{sample}/',
		cutoff = config['infercnv']['cutoff'],
		cnv_analysis_mode = config['infercnv']['cnv_analysis_mode'],
		denoise = config['infercnv']['denoise'],
		noise_logistic = config['infercnv']['noise_logistic'],
		leiden_res = config['infercnv']['leiden_res'],
	output:
		infercnv_obj = resultoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/{sample}/run.final.infercnv_obj',
		#cnv_plot = resultoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/{sample}/infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.png',
		#exprs_plot = resultoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/{sample}/infercnv.png',
	resources: 
		mem_mb = 100000
	threads: 32
	container: 
		"docker://trinityctat/infercnv"
	script:
		'cell-type-assignment/run-infercnv.R'

rule move_infercnv_figs_to_figureoutput:
	params:
		
	input:
		cnv_plot = resultoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/{sample}/infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.png',
		exprs_plot = resultoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/{sample}/infercnv.png',
	output:
		cnv_out = figureoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/{sample}/{sample}_infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.png',
		exprs_out = figureoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/{sample}/{sample}_infercnv.png',
	shell:
		'cp "{input.cnv_plot}" "{output.cnv_out}" ; cp "{input.exprs_plot}" "{output.exprs_out}"'

rule classify_malignant_cells:
	input:
		sce_tumor = dataoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/scRNASeq-infercnv-prepared-sce-{celltype}-{cohort}-ref-{normcohort}.rds',
	params:
		infercnv_output_dir = resultoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/',
		regional_cnv_calls = "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat",
		gene_cnv_calls = "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat",
		cell_groupings = "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings",
		cnv_state_var_thres = 0.92,
	output:
		sce_tumor_results_list = dataoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/scRNASeq-sce-{celltype}-{cohort}-ref-{normcohort}-results-list.rds',
		sce_tumor_infercnv_list = dataoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/scRNASeq-sce-{celltype}-{cohort}-ref-{normcohort}-infercnv-list.rds',
		cell_labels = resultoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/cell-labels-{celltype}-{cohort}-ref-{normcohort}.csv',
	resources:
		mem_mb = 10000
	threads: 16
	script:
		'cell-type-assignment/classify-malignant-cells.R'

rule plot_malignant_classifier_output:
	input:
		sce_tumor_results_list = dataoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/scRNASeq-sce-{celltype}-{cohort}-ref-{normcohort}-results-list.rds',
		sce_tumor_infercnv_list = dataoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/scRNASeq-sce-{celltype}-{cohort}-ref-{normcohort}-infercnv-list.rds',
	params:
	output:
		cell_cluster_plots = figureoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/cell-cluster-plots-{celltype}-{cohort}-ref-{normcohort}.png',
		cell_marker_exprs_plots = figureoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/cell-marker-exprs-plots-{celltype}-{cohort}-ref-{normcohort}.png',
	resources:
		mem_mb = 10000
	threads: 4
	script:
		'cell-type-assignment/plot-malignant-classifier-output.R'

rule singler_annotation:
	input:
		sce = dataoutput + 'subset-sce/{celltype}/scRNASeq-{celltype}-sce-{cohort}.rds',
		ref = config['data']['reference'] + '{celltype}-ref.rds',
	params:
		ref_label_field = lambda wildcards: singler_info.get(wildcards.celltype)
	output:
		sce = dataoutput + 'SingleR-annotation/{celltype}/scRNASeq-SingleR-annotated-{celltype}-sce-{cohort}.rds',
		dframe = resultoutput + 'SingleR-annotation/{celltype}/DFrame/{celltype}-cell-type-prediction-{cohort}.rds',
		tsv = resultoutput + 'SingleR-annotation/{celltype}/tsv/{celltype}-cell-type-prediction-scores-{cohort}.tsv',
	resources:
		mem_mb = 15000
	threads: 10
	script:
		'cell-type-assignment/singler-annotation.R'

rule visualize_singler_output:
	input:
		sce = dataoutput + 'SingleR-annotation/{celltype}/scRNASeq-SingleR-annotated-{celltype}-sce-{cohort}.rds',
		dframe = resultoutput + 'SingleR-annotation/{celltype}/DFrame/{celltype}-cell-type-prediction-{cohort}.rds',
	params:
		ncol_violin = config['singler']['ncol_for_violin_plots'],
		celltype_label_field = config['singler']['celltype_label_field'],
		individual_celltype_marker_expression_heatmap_dir = figureoutput + 'SingleR-annotation/{celltype}/celltype-marker-expression/{cohort}/',
	output:
		heatmap_assignment_score = figureoutput + 'SingleR-annotation/{celltype}/celltype-assignment-score/heatmap-{celltype}-celltype-assignment-score-{cohort}.png',
		violin_delta_distribution = figureoutput + 'SingleR-annotation/{celltype}/score-distribution/violin-{celltype}-celltype-delta-distribution-{cohort}.png',
		violin_score_distribution = figureoutput + 'SingleR-annotation/{celltype}/score-distribution/violin-{celltype}-celltype-score-distribution-{cohort}.png',
	resources:
		mem_mb = 6000
	threads: 2
	script:
		'cell-type-assignment/visualize-singler-output.R'

rule visualize_singler_output_dimred:
	input:
		sce = dataoutput + 'SingleR-annotation/{celltype}/scRNASeq-SingleR-annotated-{celltype}-sce-{cohort}.rds',
		dframe = resultoutput + 'SingleR-annotation/{celltype}/DFrame/{celltype}-cell-type-prediction-{cohort}.rds',
	params:
		celltype_label_field = config['singler']['celltype_label_field'],
		dimred_to_plot = '{dimred}',
		metadata_to_color = '{metadata}',
		sce_assay_to_plot = config['singler']['sce_assay_to_plot'],
		individual_celltype_score_dimred_dir = figureoutput + 'SingleR-annotation/{celltype}/celltype-assignment-score/individual-dimred/{cohort}/',
	output:
		dimred_cell_type = figureoutput+ 'SingleR-annotation/{celltype}/annotated-dimred/{dimred}-cell-type-and-{metadata}-{cohort}.png'
	resources: 
		mem_mb = 6000
	threads: 2
	script:
		'cell-type-assignment/visualize-singler-output-dimred.R'

rule plot_immune_cell_type_markers:
	input:
		sce = dataoutput + 'SingleR-annotation/{celltype}/scRNASeq-SingleR-annotated-{celltype}-sce-{cohort}.rds',
		dframe = resultoutput + 'SingleR-annotation/{celltype}/DFrame/{celltype}-cell-type-prediction-{cohort}.rds',
	params:
		celltype_label_field = config['singler']['celltype_label_field'],
		dimred_to_plot = config['singler']['redim_to_plot'],
		sce_assay_to_plot = config['singler']['sce_assay_to_plot'],
		marker = '{marker}',
		individual_plots_dir = figureoutput + 'SingleR-annotation/{celltype}/celltype-markers/{cohort}/individual-plots/{marker}/',
	output:
		expression_plot = figureoutput + 'SingleR-annotation/{celltype}/celltype-markers/{cohort}/{marker}.png',
	resources:
		mem_mb = 6000
	threads: 2
	script:
		'cell-type-assignment/plot-immune-cell-type-markers.R'

# rule azimuth_annotation:
# 	input:
# 		sce = dataoutput + 'process-sce/scRNASeq-filtered-sce-{cohort}.rds',
# 	params:
# 		azimuth_ref = config['azimuth']['reference'],
# 		annot_level = config['azimuth']['annotation_level'],
# 	output:
# 		seu = dataoutput + 'Azimuth-annotation/scRNASeq-Azimuth-assigned-seu-{cohort}.rds',
# 		sce = dataoutput + 'Azimuth-annotation/scRNASeq-Azimuth-assigned-sce-{cohort}.rds',
# 		assignment_score = resultoutput + 'Azimuth-annotation/cell-type-annotation-scores-{cohort}.tsv',
# 	resources:
# 		mem_mb = 10000
# 	threads: 6
# 	script:
# 		'cell-type-assignment/azimuth-annotation.R'

# rule visualize_azimuth_output:
# 	input:
# 		seu = dataoutput + 'Azimuth-annotation/scRNASeq-Azimuth-assigned-seu-{cohort}.rds',
# 		assignment_score = resultoutput + 'Azimuth-annotation/cell-type-annotation-scores-{cohort}.tsv',
# 	params:
# 		dimred_to_plot = config['azimuth']['redim_to_plot'],
# 		annot_level = config['azimuth']['annotation_level'],
# 		metadata_to_color = config['azimuth']['redim_metadata_to_color'],
# 		individual_celltype_prediction_score_UMAP_dir = figureoutput + 'Azimuth-annotation/celltype-prediction-score/individual-dimred/{cohort}/',
# 	output:
# 		umap_cell_type = figureoutput + 'Azimuth-annotation/annotated-dimred/UMAP-cell-type-and-sample-{cohort}.png',
# 		umap_prediction_score = figureoutput + 'Azimuth-annotation/celltype-prediction-score/UMAP-celltype-prediction-score-{cohort}.png',
# 		heatmap_prediction_score = figureoutput + 'Azimuth-annotation/celltype-prediction-score/heatmap-celltype-prediction-score-{cohort}.png',
# 	resources:
# 		mem_mb = 6000
# 	threads: 2
# 	script:
# 		'cell-type-assignment/visualize-azimuth-output.R'

# rule visualize_all_azimuth_output:
# 	input:
# 		sces = expand(dataoutput + 'Azimuth-annotation/scRNASeq-Azimuth-assigned-sce-{cohort}.rds', cohort = cohorts_together)
# 	params:
# 		dimred_to_plot = config['azimuth']['redim_to_plot'],
# 		annot_level = config['azimuth']['annotation_level'],
# 	output:
# 		umap_celltype = figureoutput + 'Azimuth-annotation/annotated-dimred/cohort-combined-plots/UMAP-celltype-all-cohorts.png',
# 		umap_cohort = figureoutput + 'Azimuth-annotation/annotated-dimred/cohort-combined-plots/UMAP-cohort-all-cohorts.png',
# 	resources:
# 		mem_mb = 10000
# 	threads: 2
# 	script:
# 		'cell-type-assignment/visualize-all-azimuth-output.R'

# rule plot_cell_type_markers:
# 	input: 
# 		seu = dataoutput + 'Azimuth-annotation/scRNASeq-Azimuth-assigned-seu-{cohort}.rds',
# 		sce = dataoutput + 'Azimuth-annotation/scRNASeq-Azimuth-assigned-sce-{cohort}.rds',
# 		assignment_score = resultoutput + 'Azimuth-annotation/cell-type-annotation-scores-{cohort}.tsv',
# 	params: 
# 		dimred_to_plot = config['azimuth']['redim_to_plot'],
# 		annot_level = config['azimuth']['annotation_level'],
# 		sce_assay_to_plot = config['azimuth']['sce_assay_to_plot'],
# 		marker = '{marker}',
# 		individual_plots_dir = figureoutput + 'Azimuth-annotation/celltype-markers/{cohort}/individual-plots/{marker}/', 
# 	output: 
# 		expression_plot = figureoutput + 'Azimuth-annotation/celltype-markers/{cohort}/{marker}.png',
# 	resources: 
# 		mem_mb = 3000
# 	threads: 4
# 	script: 
# 		'cell-type-assignment/plot-cell-type-markers.R'