rule subset_sce:
	input:
		sce = dataoutput + 'Azimuth-annotation/scRNASeq-Azimuth-assigned-sce-{cohort}.rds',
	params:
		annot_level = config['azimuth']['annotation_level'],
		cell_type = '{celltype}',
	output:
		sce = dataoutput + 'subset-sce/{celltype}/scRNASeq-{celltype}-sce-{cohort}.rds',
		cells_per_sample_plot = figureoutput + 'subset-sce/{celltype}/number-of-{celltype}-cells-{cohort}.png',
		cells_per_sample = resultoutput + 'subset-sce/{celltype}/number-of-{celltype}-cells-{cohort}.tsv',
	resources:
		mem_mb = 6000
	threads: 4
	script:
		'cell-selection/subset-sce.R'

rule subset_sce_further:
	input:
		sce = dataoutput + 'SingleR-annotation/{celltype}/sceRNASeq-SingleR-annotated-{celltype}-sce-{cohort}.rds',
	params:
		celltype_label_field = config['singler']['celltype_label_field'],
		subtype = '{subtype}',
	output:
		sce = dataoutput + 'subset-sce/{celltype}/{subtype}/scRNASeq-{subtype}-sce-{cohort}.rds',
		cells_per_sample_plot = figureoutput + 'subset-sce/{celltype}/{subtype}/number-of-{subtype}-cells-{cohort}.png',
		cells_per_sample = resultoutput + 'subset-sce/{celltype}/{subtype}/number-of-{subtype}-cells-{cohort}.tsv',
	resources:
		mem_mb = 4000
	threads: 4
	script:
		'cell-selection/subset-sce-further.R'

rule combine_sce:
	input:
		sces = lambda wildcards: expand(get_subsetted_sce_dir('{subtype}', dataoutput, celltype_hierarchy) + 'scRNAseq-{subtype}-sce-' + wildcards.cohort + '.rds', subtype = list(ct_zoomout_info[wildcards.ct])),
	params:
	output:
		#sce = lambda wildcards: get_combined_sce_dir(wildcards.ct, ct_zoomout_info, dataoutput, celltype_hierarchy) + wildcards.ct + '/scRNASeq-' + wildcards.ct + '-sce-' + wildcards.cohort + '.rds',
	resources:
		mem_mb = 10000
	threads: 2
	script:
		'cell-selection/combine-sce.R'