rule azimuth_annotation:
	input:
		sce = dataoutput + 'process-sce/scRNASeq-filtered-sce-{cohort}.rds',
	params:
		azimuth_ref = config['azimuth']['reference'],
		annot_level = config['azimuth']['annotation_level'],
	output:
		seu = dataoutput + 'Azimuth-annotation/scRNASeq-Azimuth-assigned-seu-{cohort}.rds',
		sce = dataoutput + 'Azimuth-annotation/scRNASeq-Azimuth-assigned-sce-{cohort}.rds',
		assignment_score = resultoutput + 'Azimuth-annotation/cell-type-annotation-scores-{cohort}.tsv',
	resources:
		mem_mb = 10000
	threads: 6
	script:
		'cell-type-assignment/azimuth-annotation.R'

rule visualize_azimuth_output:
	input:
		seu = dataoutput + 'Azimuth-annotation/scRNASeq-Azimuth-assigned-seu-{cohort}.rds',
		assignment_score = resultoutput + 'Azimuth-annotation/cell-type-annotation-scores-{cohort}.tsv',
	params:
		dimred_to_plot = config['azimuth']['redim_to_plot'],
		annot_level = config['azimuth']['annotation_level'],
		metadata_to_color = config['azimuth']['redim_metadata_to_color'],
		individual_celltype_prediction_score_UMAP_dir = figureoutput + 'Azimuth-annotation/celltype-prediction-score/individual-dimred/{cohort}/',
	output:
		umap_cell_type = figureoutput + 'Azimuth-annotation/annotated-dimred/UMAP-cell-type-and-sample-{cohort}.png',
		umap_prediction_score = figureoutput + 'Azimuth-annotation/celltype-prediction-score/UMAP-celltype-prediction-score-{cohort}.png',
		heatmap_prediction_score = figureoutput + 'Azimuth-annotation/celltype-prediction-score/heatmap-celltype-prediction-score-{cohort}.png',
	resources:
		mem_mb = 6000
	threads: 2
	script:
		'cell-type-assignment/visualize-azimuth-output.R'

rule plot_cell_type_markers:
	input: 
		seu = dataoutput + 'Azimuth-annotation/scRNASeq-Azimuth-assigned-seu-{cohort}.rds',
		sce = dataoutput + 'Azimuth-annotation/scRNASeq-Azimuth-assigned-sce-{cohort}.rds',
		assignment_score = resultoutput + 'Azimuth-annotation/cell-type-annotation-scores-{cohort}.tsv',
	params: 
		dimred_to_plot = config['azimuth']['redim_to_plot'],
		annot_level = config['azimuth']['annotation_level'],
		sce_assay_to_plot = config['azimuth']['sce_assay_to_plot'],
		marker = '{marker}',
		individual_plots_dir = figureoutput + 'Azimuth-annotation/celltype-markers/{cohort}/individual-plots/{marker}/', 
	output: 
		expression_plot = figureoutput + 'Azimuth-annotation/celltype-markers/{cohort}/{marker}.png',
	resources: 
		mem_mb = 3000
	threads: 4
	script: 
		'cell-type-assignment/plot-cell-type-markers.R'

rule run_infercnv:
	input:
		sce_normal = dataoutput + 'subset-sce/ductal/scRNASeq-{celltype}-sce-{normcohort}.rds',
		sce_tumor = dataoutput + 'subset-sce/ductal/scRNASeq-{celltype}-sce-{cohort}.rds',
		gene_order_file = config['infercnv']['gene_order_file'],
	params:
		annot_level = config['azimuth']['annotation_level'],
		infercnv_outdir = resultoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/',
	output:
		infercnv_obj = resultoutput + 'inferCNV/{celltype}/{cohort}/ref-{normcohort}/run.final.infercnv_obj',
	resources: 
		mem_mb = 25000
	threads: 16
	script:
		'cell-type-assignment/run-infercnv.R'

rule singler_annotation:
	input:
		sce = dataoutput + 'subset-sce/{celltype}/scRNASeq-{celltype}-sce-{cohort}.rds',
		ref = config['data']['reference'] + '{celltype}-ref.rds',
	params:
		ref_label_field = lambda wildcards: singler_info.get(wildcards.celltype)
	output:
		sce = dataoutput + 'SingleR-annotation/{celltype}/sceRNASeq-SingleR-annotated-{celltype}-sce-{cohort}.rds',
		dframe = resultoutput + 'SingleR-annotation/{celltype}/DFrame/{celltype}-cell-type-prediction-{cohort}.rds',
		tsv = resultoutput + 'SingleR-annotation/{celltype}/tsv/{celltype}-cell-type-prediction-scores-{cohort}.tsv',
	resources:
		mem_mb = 15000
	threads: 10
	script:
		'cell-type-assignment/singler-annotation.R'

rule visualize_singler_output:
	input:
		sce = dataoutput + 'SingleR-annotation/{celltype}/sceRNASeq-SingleR-annotated-{celltype}-sce-{cohort}.rds',
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
		sce = dataoutput + 'SingleR-annotation/{celltype}/sceRNASeq-SingleR-annotated-{celltype}-sce-{cohort}.rds',
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
		sce = dataoutput + 'SingleR-annotation/{celltype}/sceRNASeq-SingleR-annotated-{celltype}-sce-{cohort}.rds',
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