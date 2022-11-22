rule cohort_discovery_validation_grouping:
	input:
		sce = lambda wildcards: expand(get_subsetted_sce_dir(wildcards.subtype, dataoutput, celltype_hierarchy) + 'scRNAseq-' + wildcards.subtype + '-sce-{cohort}.rds', cohort = discovery_cohorts + validation_cohorts)
	params:
		cohorts_discovery = discovery_cohorts,
		cohorts_validation = validation_cohorts,
	output:
		sce_discovery = dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNAseq-{subtype}-sce-discovery.rds',
		sce_validation = dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNAseq-{subtype}-sce-validation.rds',
		scelist_discovery = dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNAseq-{subtype}-scelist-discovery.rds',
		scelist_validation = dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNAseq-{subtype}-scelist-validation.rds',
		cells_per_sample_discovery = resultoutput + 'cohort-discovery-validation-grouping/{subtype}/number-of-{subtype}-cells-discovery.tsv',
		cells_per_sample_validation = resultoutput + 'cohort-discovery-validation-grouping/{subtype}/number-of-{subtype}-cells-validation.tsv',
	resources:
		mem_mb = 15000
	threads: 10
	script:
		'data-organization/cohort-discovery-validation-grouping.R'