rule process_sce:
	input:
		sce = config['data']['sce'] + 'scRNASeq-manually-filtered-sce-{cohort}.rds',
		combined_doublet_detection = resultoutput + 'doublet-detection/DoubletFinder-combined-{cohort}.tsv',
	params:
		filtering_params = config['data_processing']['sce_filtering'],
	output: 
		sce = dataoutput + 'process-sce/scRNASeq-filtered-sce-{cohort}.rds',
		cells_per_sample_plot = figureoutput + 'process-sce/number-of-cells-{cohort}.png',
		cells_per_sample = resultoutput + 'process-sce/number-of-cells-{cohort}.tsv',
	resources:
		mem_mb = 25000
	threads: 12
	script:
		'process-data/process-sce.R'

rule process_Doublet_detection_results:
	input: 
		sce = config['data']['sce'] + 'scRNASeq-manually-filtered-sce-{cohort}.rds',
	params:
		doublet_detection_dir = config['data']['doublet_detection_results'] + '{cohort}/',
	output: 
		combined_doublet_detection = resultoutput + 'doublet-detection/DoubletFinder-combined-{cohort}.tsv',
	resources:
		mem_mb = 3000
	threads: 1
	script:
		'process-data/process-doublet-detection-results.R'