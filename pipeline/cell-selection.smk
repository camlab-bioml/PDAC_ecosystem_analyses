rule subset_sce_all_cells:
	input:
		sce = dataoutput + 'SingleR-annotation-all-cells/scRNASeq-SingleR-annotated-sce-{cohort}.rds',
	params:
		celltype_label_field = config['singler']['celltype_label_field'],
		cell_type = '{celltype}',
	output:
		sce = dataoutput + 'subset-sce/{celltype}/scRNASeq-{celltype}-sce-{cohort}.rds',
		cells_per_sample_plot = figureoutput + 'subset-sce/{celltype}/number-of-{celltype}-cells-{cohort}.png',
		cells_per_sample = resultoutput + 'subset-sce/{celltype}/number-of-{celltype}-cells-{cohort}.tsv',
	resources:
		mem_mb = 6000
	threads: 4
	script:
		'cell-selection/subset-sce-all-cells.R'

# rule subset_sce:
# 	input:
# 		sce = dataoutput + 'Azimuth-annotation/scRNASeq-Azimuth-assigned-sce-{cohort}.rds',
# 	params:
# 		annot_level = config['azimuth']['annotation_level'],
# 		cell_type = '{celltype}',
# 	output:
# 		sce = dataoutput + 'subset-sce/{celltype}/scRNASeq-{celltype}-sce-{cohort}.rds',
# 		cells_per_sample_plot = figureoutput + 'subset-sce/{celltype}/number-of-{celltype}-cells-{cohort}.png',
# 		cells_per_sample = resultoutput + 'subset-sce/{celltype}/number-of-{celltype}-cells-{cohort}.tsv',
# 	resources:
# 		mem_mb = 6000
# 	threads: 4
# 	script:
# 		'cell-selection/subset-sce.R'

rule subset_sce_further:
	input:
		sce = dataoutput + 'SingleR-annotation/{celltype}/scRNASeq-SingleR-annotated-{celltype}-sce-{cohort}.rds',
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
		sces = lambda wildcards: expand(dataoutput + 'subset-sce/{subtype}/scRNASeq-{subtype}-sce-' + wildcards.cohort + '.rds', subtype = list(ct_zoomout_info[wildcards.ct])),
	params:
		celltype_label_field = config['singler']['celltype_label_field'],
		cell_type = '{ct}',
	output:
		sce = dataoutput + 'merged-sce/{ct}/scRNASeq-{ct}-sce-{cohort}.rds',
		cells_per_sample_plot = figureoutput + 'merged-sce/{ct}/number-of-{ct}-cells-{cohort}.png',
		cells_per_sample = resultoutput + 'merged-sce/{ct}/number-of-{ct}-cells-{cohort}.tsv',
	resources:
		mem_mb = 15000
	threads: 4
	script:
		'cell-selection/combine-sce.R'