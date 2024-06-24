rule cohort_discovery_validation_grouping_and_sample_filtering:
	input:
		sce = lambda wildcards: expand(get_sce_dir(wildcards.subtype, dataoutput, celltypes, ct_zoomout, celltype_hierarchy) + 'scRNASeq-' + wildcards.subtype + '-sce-{cohort}.rds', cohort = discovery_cohorts + validation_cohorts)
	params:
		cohorts_discovery = discovery_cohorts,
		cohorts_validation = validation_cohorts,
		sample_cell_count_thres = config['liger']['sample_cell_count_thres'],
	output:
		sce_discovery = dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNASeq-{subtype}-sce-discovery.rds',
		sce_validation = dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNASeq-{subtype}-sce-validation.rds',
		scelist_discovery = dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNASeq-{subtype}-scelist-discovery.rds',
		scelist_validation = dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNASeq-{subtype}-scelist-validation.rds',
		cells_per_sample_discovery = resultoutput + 'cohort-discovery-validation-grouping/{subtype}/number-of-{subtype}-cells-discovery.tsv',
		cells_per_sample_validation = resultoutput + 'cohort-discovery-validation-grouping/{subtype}/number-of-{subtype}-cells-validation.tsv',
	resources:
		mem_mb = 50000
	threads: 20
	script:
		'data-organization/cohort-discovery-validation-grouping.R'

rule cohort_examination:
	input:
		sces = lambda wildcards: expand(dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNASeq-{subtype}-sce-' + wildcards.group + '.rds', subtype = ct_extract),
	params:
		cohorts_discovery = discovery_cohorts,
		cohorts_validation = validation_cohorts,
		celltypes = ct_extract,
	output:
		sce = dataoutput + 'cohort-discovery-validation-grouping/cohort-examination/sce-prepared-and-sampled-{group}.rds',
	resources:
		mem_mb = 100000
	threads: 24
	script:
		'data-organization/cohort-examination.R'

rule draw_cohort_examination:
	input:
		sce = dataoutput + 'cohort-discovery-validation-grouping/cohort-examination/sce-prepared-and-sampled-{group}.rds',
	params:
		cohorts_discovery = discovery_cohorts,
		cohorts_validation = validation_cohorts,
		rle_style = "minimal",
	output:
		rle_plot_counts = figureoutput + 'cohort-discovery-validation-grouping/cohort-examination/RLE-plot-counts-{group}.png',
		rle_plot_logcounts = figureoutput + 'cohort-discovery-validation-grouping/cohort-examination/RLE-plot-logcounts-{group}.png',
		rle_plot_seuratNormData = figureoutput + 'cohort-discovery-validation-grouping/cohort-examination/RLE-plot-seuratNormData-{group}.png',
		logcounts_mean_plot = figureoutput + 'cohort-discovery-validation-grouping/cohort-examination/logcounts-mean-plot-{group}.png',
		logcounts_95_perct_plot = figureoutput + 'cohort-discovery-validation-grouping/cohort-examination/logcounts-95-perct-plot-{group}.png',
	resources:
		mem_mb = 4000
	threads: 1
	script:
		'data-organization/draw-cohort-examination.R'

rule figure_1_preparation:
	input:
		sces = lambda wildcards: expand(dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNASeq-{subtype}-sce-' + wildcards.group + '.rds', subtype = ct_extract),
	params:
		rpath = config['RPATH'],
		cohorts_discovery = discovery_cohorts,
		cohorts_validation = validation_cohorts,
		celltypes = ct_extract,
		dim_red_plot = config['figure1']['dim_red_plot'],
	output:
		metadata = resultoutput + 'cohort-discovery-validation-grouping/figure-1/metadata-{group}.tsv',
		sce = dataoutput + 'cohort-discovery-validation-grouping/figure-1/sce-prepared-{group}.rds',
		dimred = resultoutput + 'cohort-discovery-validation-grouping/figure-1/dimred-{group}.tsv',
	resources:
		mem_mb = 120000
	threads: 32
	script:
		'data-organization/figure-1-preparation.R'

rule draw_figure_1:
	input:
		cell_type_rename = config['figure1']['cell_type_rename_csv'],
		schematic = config['figure1']['schematic_plot'],
		metadata_dis = resultoutput + 'cohort-discovery-validation-grouping/figure-1/metadata-discovery.tsv',
		dimred_dis = resultoutput + 'cohort-discovery-validation-grouping/figure-1/dimred-discovery.tsv',
		metadata_val = resultoutput + 'cohort-discovery-validation-grouping/figure-1/metadata-validation.tsv',
		dimred_val = resultoutput + 'cohort-discovery-validation-grouping/figure-1/dimred-validation.tsv',
		sce_dis = dataoutput + 'cohort-discovery-validation-grouping/figure-1/sce-prepared-discovery.rds',
		sce_val = dataoutput + 'cohort-discovery-validation-grouping/figure-1/sce-prepared-validation.rds',
	params:
		cell_type_pallete = config['figure1']['cell_type_pallete_to_use'],
		cohorts_discovery = discovery_cohorts,
		cohorts_validation = validation_cohorts,
		metadata_plot_width = config['figure1']['metadata_plot_width'],
  		metadata_plot_height = config['figure1']['metadata_plot_height'],
  		umap_plot_width = config['figure1']['umap_plot_width'],
  		umap_plot_height = config['figure1']['umap_plot_height'],
		stacked_bar_plot_width = config['figure1']['stacked_bar_plot_width'],
		stacked_bar_plot_height = config['figure1']['stacked_bar_plot_height'],
		marker_dot_plot_width = config['figure1']['marker_dot_plot_width'],
		marker_dot_plot_height = config['figure1']['marker_dot_plot_height'],
  		figure1_width = config['figure1']['figure1_width'],
  		figure1_height = config['figure1']['figure1_height'],
	output:
		cohort_pal = figureoutput + 'cohort-palette.rds',
		celltype_pal = figureoutput + 'celltype-palette.rds',
		figure1_a = figureoutput + 'cohort-discovery-validation-grouping/figure-1/figure-1-A.png',
		figure1_c = figureoutput + 'cohort-discovery-validation-grouping/figure-1/figure-1-C.png',
		figure1_d = figureoutput + 'cohort-discovery-validation-grouping/figure-1/figure-1-D.png',
		figure1_e = figureoutput + 'cohort-discovery-validation-grouping/figure-1/figure-1-E.png',
		figure1_f = figureoutput + 'cohort-discovery-validation-grouping/figure-1/figure-1-F.png',
		figure1_png = figureoutput + 'cohort-discovery-validation-grouping/figure-1/figure-1.png',
		figure1_pdf = figureoutput + 'cohort-discovery-validation-grouping/figure-1/figure-1.pdf',
	resources:
		mem_mb = 80000
	threads: 32
	script:
		'data-organization/draw-figure-1.R'