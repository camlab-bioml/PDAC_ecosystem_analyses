#shell.executable("/bin/bash")
#shell.prefix("source ~/.bashrc; ")

### [ SET UP ] #####
import pandas as pd
import numpy as np 
import yaml
from itertools import chain
import os

container: "pdactme_ver_3.sif"
configfile: 'config/config.yml'

# Functions 
def get_mem_mb(attempt, base_mem):
    return base_mem + (attempt * 2000)

def get_sce_dir(celltype, dataoutput_dir, ct_normal, ct_broad, celltype_hierarchy):	
	if celltype in ct_normal:
		return dataoutput_dir + 'subset-sce/' + celltype + '/'

	if celltype in ct_broad:
		return dataoutput_dir + 'merged-sce/' + celltype + '/'

	for ct in celltype_hierarchy.keys():
		if celltype in celltype_hierarchy.get(ct):
			return dataoutput_dir + 'subset-sce/' + ct + '/' + celltype + '/'
		else:
			return None

def get_liger_param(parameterlist, celltype, param):
	ct_parameterlist = parameterlist[parameterlist['celltype'] == celltype]
	return list(ct_parameterlist[param])[0]

def get_subsetted_sce_dir(celltype, dataoutput_dir, celltype_hierarchy):
	root_dir = dataoutput_dir + 'subset-sce/'
	
	if celltype in celltype_hierarchy.keys():
		return root_dir + celltype + '/'

	for ct in celltype_hierarchy.keys():
		if celltype in celltype_hierarchy.get(ct):
			return root_dir + ct + '/' + celltype + '/'
		else:
			return None

def get_combined_sce_dir(celltype, celltype_info, dataoutput_dir, celltype_hierarchy):
	root_dir = dataoutput_dir + 'subset-sce/'

	ct_child = list(celltype_info[celltype])[0]
	
	if ct_child in celltype_hierarchy.keys():
		return root_dir 

	for ct in celltype_hierarchy.keys():
		if ct_child in celltype_hierarchy.get(ct):
			return root_dir + ct + '/'
		else:
			return None

# get cohorts for annotation and subsetting
cohorts_csv = pd.read_csv(config['cohorts'], header = 0)
cohorts = list(cohorts_csv['cohort'])

# get cohorts for discovery and validation grouping
if (config['liger']['discovery_cohort_list'] is not None) & os.path.isfile(config['liger']['discovery_cohort_list']):
	discovery_cohorts = pd.read_csv(config['liger']['discovery_cohort_list'], header = 0)
	discovery_cohorts = list(discovery_cohorts['cohort'])
	validation_cohorts = pd.read_csv(config['liger']['validation_cohort_list'], header = 0)
	validation_cohorts = list(validation_cohorts['cohort'])

# This determines the output folder - the version can be set in the config file
output = 'output/' + config['version'] + '/'
figureoutput = output + 'figures/'
dataoutput = output + 'data/'
resultoutput = output + 'results/'
groupedfilesoutput = output + 'groupedfiles/'



# Get list of markers
markers = yaml.safe_load(open(config['pdac_cell_type_markers']))
unique_markers = []
for i in markers:
    for j in markers[i]:
        unique_markers.append(markers[i][j])

unique_markers = list(set(list(chain.from_iterable(unique_markers))))

immune_markers = yaml.safe_load(open(config['immune_cell_type_markers']))
unique_immune_markers = []
for i in immune_markers:
	for j in immune_markers[i]:
		for h in immune_markers[i][j]:
			unique_immune_markers.append(immune_markers[i][j][h])	

unique_immune_markers = list(set(list(chain.from_iterable(unique_immune_markers))))

celltype_hierarchy = dict()

### [ LIST OUTPUTS TO BE CREATED ] #####
# grouped outputs
groupedfiles = {
	# 'combined_sce': expand(groupedfilesoutput + 'combined-sce/combined-sce-{cohort}.rds', cohort = cohorts),
	# 'combined_sce_metadata': expand(groupedfilesoutput + 'combined-sce/combined-sce-metadata-{cohort}.tsv', cohort = cohorts),
	# 'combined_sce_metadata_plot': expand(figureoutput + 'combined-sce/combined-sce-metadata-{cohort}.png', cohort = cohorts),
	# 'combined_sce_metadata_plot_annotated': expand(figureoutput + 'combined-sce/combined-sce-metadata-annotated-{cohort}.png', cohort = cohorts),
	# 'combined_sce_metadata_plot_annotated_with_celltype': expand(figureoutput + 'combined-sce/combined-sce-metadata-annotated-with-celltype-{cohort}.png', cohort = cohorts),
	# 'combined_sce_metadata_plot_annotated_with_celltype_and_sample': expand(figureoutput + 'combined-sce/combined-sce-metadata-annotated-with-celltype-and-sample-{cohort}.png', cohort = cohorts),
	# 'combined_sce_metadata_plot_annotated_with_celltype_and_sample_and_immune': expand(figureoutput + 'combined-sce/combined-sce-metadata-annotated-with-celltype-and-sample-and-immune-{cohort}.png', cohort = cohorts),
	# 'combined_sce_metadata_plot_annotated_with_celltype_and_sample_and_immune_and_cellcycle': expand(figureoutput + 'combined-sce/combined-sce-metadata-annotated-with-celltype-and-sample-and-immune-and-cellcycle-{cohort}.png', cohort = cohorts),
	# 'combined_sce_metadata_plot_annotated_with_celltype_and_sample_and_immune_and_cellcycle_and_doublet': expand(figureoutput + 'combined-sce/combined-sce-metadata-annotated-with-celltype-and-sample-and-immune-and-cellcycle-and-doublet-{cohort}.png', cohort = cohorts),
	# 'combined_sce_metadata_plot_annotated_with_celltype_and_sample_and_immune_and_cellcycle_and_doublet_and_pseudotime': expand(figureoutput + 'combined-sce/combined-sce-metadata-annotated-with-celltype-and-sample-and-immune-and-cellcycle-and-doublet-and-pseudotime-{cohort}.png', cohort = cohorts),
	# 'combined_sce_metadata_plot_annotated_with_celltype_and_sample_and_immune_and_cellcycle_and_doublet_and_pseudotime_and_Azimuth': expand(figureoutput + 'combined-sce/combined-sce-metadata-annotated-with-celltype-and-sample-and-immune-and-cellcycle-and-doublet-and-pseudotime-and-Azimuth-{cohort}.png', cohort = cohorts),
}
# outputs
results = {
	'combined_doublet_detection': expand(resultoutput + 'doublet-detection/DoubletFinder-combined-{cohort}.tsv', cohort = cohorts),

    	'scRNA_filtered_sce': expand(dataoutput + 'process-sce/scRNASeq-filtered-sce-{cohort}.rds', cohort = cohorts),
	'cells_per_sample': expand(resultoutput + 'process-sce/number-of-cells-{cohort}.tsv', cohort = cohorts),
    	'cells_per_sample_plot': expand(figureoutput + 'process-sce/number-of-cells-{cohort}.png', cohort = cohorts),

	#'scRNA_Azimuth_assigned_seu': expand(dataoutput + 'Azimuth-annotation/scRNASeq-Azimuth-assigned-seu-{cohort}.rds', cohort = cohorts),
    	#'scRNA_Azimuth_assigned_sce': expand(dataoutput + 'Azimuth-annotation/scRNASeq-Azimuth-assigned-sce-{cohort}.rds', cohort = cohorts),
	#'cell_type_assignment_scores': expand(resultoutput + 'Azimuth-annotation/cell-type-annotation-scores-{cohort}.tsv', cohort = cohorts),

	'scRNA_SingleR_annotated_all_cells_sce': expand(dataoutput + 'SingleR-annotation-all-cells/scRNASeq-SingleR-annotated-sce-{cohort}.rds', cohort = cohorts),
	'scRNA_SingleR_all_cell_prediction_dframe': expand(resultoutput + 'SingleR-annotation-all-cells/DFrame/cell-type-prediction-{cohort}.rds', cohort = cohorts),

	'SingleR_assignment_score_heatmap_all_cells': expand(figureoutput + 'SingleR-annotation-all-cells/celltype-assignment-score/heatmap-celltype-assignment-score-{cohort}.png', cohort = cohorts),
	'SingleR_delta_distribution_violin_all_cells': expand(figureoutput + 'SingleR-annotation-all-cells/score-distribution/violin-celltype-delta-distribution-{cohort}.png', cohort = cohorts),
	'SingleR_score_distribution_violin_all_cells': expand(figureoutput + 'SingleR-annotation-all-cells/score-distribution/violin-celltype-score-distribution-{cohort}.png', cohort = cohorts),

	'SingleR_cell_type_and_sample_dimred_plot_all_cells': expand(figureoutput+ 'SingleR-annotation-all-cells/annotated-dimred/{dimred}-cell-type-and-{metadata}-{cohort}.png', cohort = cohorts, dimred = config['singler']['redim_to_plot'], metadata = config['singler']['redim_metadata_to_color']),

	#'cell_type_and_sample_UMAP_plot': expand(figureoutput + 'Azimuth-annotation/annotated-dimred/UMAP-cell-type-and-sample-{cohort}.png', cohort = cohorts),
	#'cell_type_prediction_score_UMAP_plot': expand(figureoutput + 'Azimuth-annotation/celltype-prediction-score/UMAP-celltype-prediction-score-{cohort}.png', cohort = cohorts),
	#'cell_type_prediction_score_heatmap': expand(figureoutput + 'Azimuth-annotation/celltype-prediction-score/heatmap-celltype-prediction-score-{cohort}.png', cohort = cohorts),
}

# if (config['visualize_these_cohorts_together'] is not None) & os.path.isfile(config['visualize_these_cohorts_together']): 
# 	# Visualize cell type labels predicted by Azimuth for all cohorts in one plot
# 	cohorts_together = pd.read_csv(config['visualize_these_cohorts_together'], header = 0)
# 	cohorts_together = list(cohorts_together['cohort'])

# 	colour_field = ['cohort', 'celltype']

# 	results['combined_UMAP_plots'] = expand(figureoutput + 'Azimuth-annotation/annotated-dimred/cohort-combined-plots/UMAP-{field}-all-cohorts.png', field = colour_field) 

# if unique_markers:
# 	results['cell_type_marker_plot-all-cells'] = expand(figureoutput + 'SingleR-annotation-all-cells/celltype-markers/{cohort}/{marker}.png', cohort = cohorts, marker = unique_markers)

if (config['celltypes_csv'] is not None) & os.path.isfile(config['celltypes_csv']):
	# Create a list of all cell types of interest
	celltypes = pd.read_csv(config['celltypes_csv'], header = 0) 
	celltypes = list(celltypes['celltype'])
	
	celltype_hierarchy = {k: [] for k in celltypes}

	results['scRNA_subsetted_sce'] = expand(dataoutput + 'subset-sce/{celltype}/scRNASeq-{celltype}-sce-{cohort}.rds', cohort = cohorts, celltype = celltypes)
	results['subsetted_cells_per_sample_plot'] = expand(figureoutput + 'subset-sce/{celltype}/number-of-{celltype}-cells-{cohort}.png', cohort = cohorts, celltype = celltypes)
	results['subsetted_cells_per_sample'] = expand(resultoutput + 'subset-sce/{celltype}/number-of-{celltype}-cells-{cohort}.tsv', cohort = cohorts, celltype = celltypes)

if (config['singler']['celltypes_focused'] is not None) & os.path.isfile(config['singler']['celltypes_focused']):
	# Create a list of cell types for subtype annotation
	singler = pd.read_csv(config['singler']['celltypes_focused'], header = 0)
	singler_info = dict(zip(singler['celltype'], singler['ref_label_field']))
	ct_zoomin = singler_info.keys()

	results['scRNA_SingleR_annotated_sce'] = expand(dataoutput + 'SingleR-annotation/{celltype}/scRNASeq-SingleR-annotated-{celltype}-sce-{cohort}.rds', cohort = cohorts, celltype = ct_zoomin)
	results['zoom_in_cell_type_prediction_dframe'] = expand(resultoutput + 'SingleR-annotation/{celltype}/DFrame/{celltype}-cell-type-prediction-{cohort}.rds', cohort = cohorts, celltype = ct_zoomin)
	results['zoom_in_cell_type_prediction_tsv'] = expand(resultoutput + 'SingleR-annotation/{celltype}/tsv/{celltype}-cell-type-prediction-scores-{cohort}.tsv', cohort = cohorts, celltype = ct_zoomin)

	results['SingleR_assignment_score_heatmap'] = expand(figureoutput + 'SingleR-annotation/{celltype}/celltype-assignment-score/heatmap-{celltype}-celltype-assignment-score-{cohort}.png', cohort = cohorts, celltype = ct_zoomin)
	results['SingleR_delta_distribution_violin'] = expand(figureoutput + 'SingleR-annotation/{celltype}/score-distribution/violin-{celltype}-celltype-delta-distribution-{cohort}.png', cohort = cohorts, celltype = ct_zoomin)
	results['SingleR_score_distribution_violin'] = expand(figureoutput + 'SingleR-annotation/{celltype}/score-distribution/violin-{celltype}-celltype-score-distribution-{cohort}.png', cohort = cohorts, celltype = ct_zoomin)

	if (config['singler']['redim_to_plot'] is not None) & (config['singler']['redim_metadata_to_color'] is not None):
		results['SingleR_cell_type_and_sample_dimred_plot'] = expand(figureoutput+ 'SingleR-annotation/{celltype}/annotated-dimred/{dimred}-cell-type-and-{metadata}-{cohort}.png', cohort = cohorts, celltype = ct_zoomin, dimred = config['singler']['redim_to_plot'], metadata = config['singler']['redim_metadata_to_color']),

	for ct in ct_zoomin:
		if os.path.isfile(config['celltype_annot_resource_dir'] + ct + '_cell_types.csv'):
			subtypes = pd.read_csv(config['celltype_annot_resource_dir'] + ct + '_cell_types.csv', header = 0)
			subtypes = list(subtypes['celltype'])

			celltype_hierarchy[ct] = subtypes

			results['scRNA_subsetted_sce'].extend(expand(dataoutput + 'subset-sce/' + ct + '/{subtype}/scRNASeq-{subtype}-sce-{cohort}.rds', cohort = cohorts, subtype = subtypes))
			results['subsetted_cells_per_sample_plot'].extend(expand(figureoutput + 'subset-sce/' + ct + '/{subtype}/number-of-{subtype}-cells-{cohort}.png', cohort = cohorts, subtype = subtypes))
			results['subsetted_cells_per_sample'].extend(expand(resultoutput + 'subset-sce/' + ct + '/{subtype}/number-of-{subtype}-cells-{cohort}.tsv', cohort = cohorts, subtype = subtypes))

# 		if ct == "immune":
# 			results['immune_cell_type_marker_plot'] = expand(figureoutput + 'SingleR-annotation/{celltype}/celltype-markers/{cohort}/{marker}.png', cohort = cohorts, celltype = ct, marker = unique_immune_markers),

if os.path.isfile(config['celltype_annot_resource_dir'] + 'cell_types_to_combine.csv'):
	ct_zoomout_info = pd.read_csv(config['celltype_annot_resource_dir'] + 'cell_types_to_combine.csv', header = 0)
	ct_zoomout = list(ct_zoomout_info.keys())
	results['scRNA_combined_sce'] = []

	for ct in ct_zoomout:
		results['scRNA_combined_sce'].extend(expand(dataoutput + 'merged-sce/' + ct + '/scRNASeq-' + ct + '-sce-{cohort}.rds', cohort = cohorts))

if config['infercnv']['run_infercnv'] & (config['infercnv']['gene_order_file'] is not None) & os.path.isfile(config['infercnv']['gene_order_file']):
	gene_order_file_path = config['infercnv']['gene_order_file']
	ct_infercnv = pd.read_csv(config['infercnv']['celltypes_for_infercnv'], header = 0)
	ct_infercnv = list(ct_infercnv['celltype'])
	cohorts_infercnv = pd.read_csv(config['infercnv']['run_infercnv_cohorts'], header = 0)
	cohorts_infercnv = list(cohorts_infercnv['cohort'])

	#results['inferCNV_sce'] = expand(dataoutput + 'inferCNV/{celltype}/{cohort}/ref-' + config['infercnv']['ref_cohort'] + '/scRNASeq-infercnv-prepared-sce-{celltype}-{cohort}-ref-' + config['infercnv']['ref_cohort'] + '.rds', celltype = ct_infercnv, cohort = cohorts_infercnv)
	results['inferCNV_objects'] = []
	results['inferCNV_cnv_plots'] = []
	results['inferCNV_exprs_plots'] = []
	results['sce_tumor_results_list'] = []
	results['sce_tumor_infercnv_list'] = []
	results['cell_labels'] = []
	results['cell_cluster_plots'] = []
	results['cell_marker_exprs_plots'] = []

	for celltype in ct_infercnv:
		for cohort in cohorts_infercnv:
			if os.path.isfile(resultoutput + 'inferCNV/' + celltype + '/' + cohort + '/ref-' + config['infercnv']['ref_cohort'] + '/sample-list.csv'):
				samples = pd.read_csv(resultoutput + 'inferCNV/' + celltype + '/' + cohort + '/ref-' + config['infercnv']['ref_cohort'] + '/sample-list.csv', header = 0)
				samples = list(samples['sample'])

				#results['inferCNV_objects'].extend(expand(resultoutput + 'inferCNV/' + celltype + '/' + cohort + '/ref-' + config['infercnv']['ref_cohort'] + '/{sample}/run.final.infercnv_obj', sample = samples))

				#results['inferCNV_cnv_plots'].extend(expand(figureoutput + 'inferCNV/' + celltype + '/' + cohort + '/ref-' + config['infercnv']['ref_cohort'] + '/{sample}/{sample}_infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.png', sample = samples))
				#results['inferCNV_exprs_plots'].extend(expand(figureoutput + 'inferCNV/' + celltype + '/' + cohort + '/ref-' + config['infercnv']['ref_cohort'] + '/{sample}/{sample}_infercnv.png', sample = samples))

				#results['sce_tumor_results_list'].extend([dataoutput + 'inferCNV/' + celltype + '/' + cohort + '/ref-' + config['infercnv']['ref_cohort'] + '/scRNASeq-sce-' + celltype + '-' + cohort + '-ref-' + config['infercnv']['ref_cohort'] + '-results-list.rds'])
				#results['sce_tumor_infercnv_list'].extend([dataoutput + 'inferCNV/' + celltype + '/' + cohort + '/ref-' + config['infercnv']['ref_cohort'] + '/scRNASeq-sce-' + celltype + '-' + cohort + '-ref-' + config['infercnv']['ref_cohort'] + '-infercnv-list.rds'])
				#results['cell_labels'].extend([resultoutput + 'inferCNV/' + celltype + '/' + cohort + '/ref-' + config['infercnv']['ref_cohort'] + '/cell-labels-' + celltype + '-' + cohort + '-ref-' + config['infercnv']['ref_cohort'] + '.csv'])

				#results['cell_cluster_plots'].extend([figureoutput + 'inferCNV/' + celltype + '/' + cohort + '/ref-' + config['infercnv']['ref_cohort'] + '/cell-cluster-plots-' + celltype + '-' + cohort + '-ref-' + config['infercnv']['ref_cohort'] + '.png'])
				#results['cell_marker_exprs_plots'].extend([figureoutput + 'inferCNV/' + celltype + '/' + cohort + '/ref-' + config['infercnv']['ref_cohort'] + '/cell-marker-exprs-plots-' + celltype + '-' + cohort + '-ref-' + config['infercnv']['ref_cohort'] + '.png'])

if (config['liger']['celltypes_for_signature'] is not None) & os.path.isfile(config['liger']['celltypes_for_signature']):
	# Create a list of cell types for extracting signatures
	ct_signature = pd.read_csv(config['liger']['celltypes_for_signature'], header = 0)
	ct_signature = list(ct_signature['celltype'])
	
	if (config['liger']['discovery_cohort_list'] is not None) & os.path.isfile(config['liger']['discovery_cohort_list']):
		discovery_cohorts = pd.read_csv(config['liger']['discovery_cohort_list'], header = 0)
		discovery_cohorts = list(discovery_cohorts['cohort'])
		validation_cohorts = pd.read_csv(config['liger']['validation_cohort_list'], header = 0)
		validation_cohorts = list(validation_cohorts['cohort'])

		results['scRNA_discovery_sce_combined'] = expand(dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNASeq-{subtype}-sce-discovery.rds', subtype = ct_signature)
		results['scRNA_validation_sce_combined'] = expand(dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNASeq-{subtype}-sce-validation.rds', subtype = ct_signature)
		results['scRNA_discovery_sce_list'] = expand(dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNASeq-{subtype}-scelist-discovery.rds', subtype = ct_signature)
		results['scRNA_validation_sce_list'] = expand(dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNASeq-{subtype}-scelist-validation.rds', subtype = ct_signature)
		
		results['discovery_cells_per_sample'] = expand(resultoutput + 'cohort-discovery-validation-grouping/{subtype}/number-of-{subtype}-cells-discovery.tsv', subtype = ct_signature)
		results['validation_cells_per_sample'] = expand(resultoutput + 'cohort-discovery-validation-grouping/{subtype}/number-of-{subtype}-cells-validation.tsv', subtype = ct_signature)

		conditions = ['discovery', 'validation']

# 		if (config['liger']['seedlist'] is not None) & os.path.isfile(config['liger']['seedlist']):
# 			seeds = pd.read_csv(config['liger']['seedlist'], header = 0)
# 			seeds = list(seeds['seed'])

# 			if (config['liger']['parameter_range'] is not None) & os.path.isfile(config['liger']['parameter_range']):
# 				parameter_range = pd.read_csv(config['liger']['parameter_range'], header = 0)
# 				lambda_range = list(parameter_range['lambda'])
# 				k_range = list(parameter_range['k'])
# 				lambdalist = list(range(lambda_range[0], lambda_range[1]))
# 				klist = list(range(k_range[0], k_range[1]))
# 			else:
# 				lambdalist = np.arange(1, 5, 1).tolist()
# 				klist = np.arange(5, 20, 5).tolist()

# 			results['LIGER_parameter_sweep_metrics'] = expand(resultoutput + 'LIGER/parameter-sweep/metrics/{subtype}/{condition}/seed-{seed}/{subtype}-{condition}-seed-{seed}-lambda-{Lambda}-k-{k}.tsv', subtype = ct_signature, condition = conditions, seed = seeds, Lambda = lambdalist, k = klist)
# 			results['LIGER_parameter_sweep_factor_contributions'] = expand(resultoutput + 'LIGER/parameter-sweep/factor-contributions/{subtype}/{condition}/seed-{seed}/k-{k}/{subtype}-{condition}-seed-{seed}-lambda-{Lambda}-k-{k}.tsv', subtype = ct_signature, condition = conditions, seed = seeds, Lambda = lambdalist, k = klist)

# 			results['LIGER_parameter_sweep_metrics_combined'] = expand(resultoutput + 'LIGER/parameter-sweep/metrics/{subtype}/{condition}/{subtype}-{condition}-parameter-sweep-metrics.tsv', subtype = ct_signature, condition = conditions)
# 			results['LIGER_parameter_sweep_metrics_plot_k'] = expand(figureoutput + 'LIGER/parameter-sweep/metrics/{subtype}/{condition}/{subtype}-{condition}-k-sweep.png', subtype = ct_signature, condition = conditions)
# 			results['LIGER_parameter_sweep_metrics_plot_lambda'] = expand(figureoutput + 'LIGER/parameter-sweep/metrics/{subtype}/{condition}/{subtype}-{condition}-lambda-sweep.png', subtype = ct_signature, condition = conditions)

		if (config['liger']['parameterlist'] is not None) & os.path.isfile(config['liger']['parameterlist']):
			parameterlist = pd.read_csv(config['liger']['parameterlist'], header = 0)
			ct_extract = list(parameterlist['celltype'])
			
			if (config['figure1']['make_figure1']): 
				results['figure_1_metadata'] = expand(resultoutput + 'cohort-discovery-validation-grouping/figure-1/metadata-{group}.tsv', group = conditions)
				results['figure_1_sce'] = expand(dataoutput + 'cohort-discovery-validation-grouping/figure-1/sce-prepared-{group}.rds', group = conditions)
				results['figure_1_dimred'] = expand(resultoutput + 'cohort-discovery-validation-grouping/figure-1/dimred-{group}.tsv', group = conditions)
				results['figure_1_plot'] = expand(figureoutput + 'cohort-discovery-validation-grouping/figure-1/figure-1.png')

				results['cohort_examination_sce'] = expand(dataoutput + 'cohort-discovery-validation-grouping/cohort-examination/sce-prepared-and-sampled-{group}.rds', group = conditions)
				#results['cohort_examination_counts'] = expand(resultoutput + 'cohort-discovery-validation-grouping/cohort-examination/counts-{group}.tsv', group = conditions)
				results['cohort_examination_RLE_plot_logcounts'] = expand(figureoutput + 'cohort-discovery-validation-grouping/cohort-examination/RLE-plot-logcounts-{group}.png', group = conditions)
				results['cohort_examination_RLE_plot_seuratNormData'] = expand(figureoutput + 'cohort-discovery-validation-grouping/cohort-examination/RLE-plot-seuratNormData-{group}.png', group = conditions)
			
			results['LIGER_variable_gene_lists'] = expand(resultoutput + 'LIGER/signature-extraction/variable-gene-list/{subtype}/{subtype}-variable-genes-{condition}.tsv', subtype = ct_extract, condition = conditions)
			results['LIGER_common_variable_gene_lists'] = expand(resultoutput + 'LIGER/signature-extraction/variable-gene-list/{subtype}/{subtype}-variable-genes-common.tsv', subtype = ct_extract)

			results['LIGER_signature_extraction_results'] = expand(resultoutput + 'LIGER/signature-extraction/LIGER-object/{subtype}/{subtype}-liger-{condition}.rds', subtype = ct_extract, condition = conditions)

			results['LIGER_signature_loading_plot'] = expand(figureoutput + 'LIGER/signature-extraction/{subtype}/{condition}/{subtype}-{condition}-signature-loading.png', subtype = ct_extract, condition = conditions)
			results['LIGER_gene_loading_plot'] = expand(figureoutput + 'LIGER/signature-extraction/{subtype}/{condition}/{subtype}-{condition}-gene-loading.png', subtype = ct_extract, condition = conditions)
			results['LIGER_dataset_cluster_plot'] = expand(figureoutput + 'LIGER/signature-extraction/{subtype}/{condition}/{subtype}-{condition}-dataset-cluster.png', subtype = ct_extract, condition = conditions)

			results['LIGER_signature_loading_matrix'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-{condition}.tsv', subtype = ct_extract, condition = conditions)
			results['LIGER_gene_loading_matrix_combined'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading.tsv', subtype = ct_extract)

			results['LIGER_gene_loading_corr'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-corr.tsv', subtype = ct_extract)
			results['LIGER_gene_loading_dist'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-dist.tsv', subtype = ct_extract)
			results['LIGER_validated_sig_df'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-validated-signatures.tsv', subtype = ct_extract)
			results['LIGER_gene_loading_corr_plot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-corr.png', subtype = ct_extract)
			results['LIGER_gene_loading_dist_plot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-dist.png', subtype = ct_extract)

			results['LIGER_validated_sig_df_alternative'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-validated-signatures-hungarian-method.tsv', subtype = ct_extract)
			results['LIGER_gene_loading_corr_plot_alternative'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-corr-hungarian-method.png', subtype = ct_extract)
			results['LIGER_gene_loading_dist_plot_alternative'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-dist-hungarian-method.png', subtype = ct_extract)
			results['LIGER_gene_loading_corr_and_dist_plot_alternative'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-corr-and-dist-hungarian-method.png', subtype = ct_extract)

			results['LIGER_simulated_signature_distances'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-simulated-distance-list.tsv', subtype = ct_extract)
			results['LIGER_simulated_signature_correlation'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-simulated-correlation-list.tsv', subtype = ct_extract)
			results['LIGER_simulated_signature_distances_plot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-simulated-distance.pdf', subtype = ct_extract)
			results['LIGER_simulated_signature_correlation_plot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-validation/{subtype}-gene-loading-simulated-correlation.pdf', subtype = ct_extract)

			results['LIGER_validated_signature_loading_matrix'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-validated.tsv', subtype = ct_extract)
			results['LIGER_validated_gene_loading_matrix'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-validated.tsv', subtype = ct_extract)
			
			results['LIGER_validated_signature_loading_correlation'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-collapse/{subtype}-signature-loading-correlation-validated.tsv', subtype = ct_extract)
			results['LIGER_validated_gene_loading_correlation'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-collapse/{subtype}-gene-loading-correlation-validated.tsv', subtype = ct_extract)

			results['LIGER_validated_signature_loading_correlation_plot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-collapse/{subtype}-signature-loading-correlation-validated.png', subtype = ct_extract)
			results['LIGER_validated_gene_loading_correlation_plot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-collapse/{subtype}-gene-loading-correlation-validated.png', subtype = ct_extract)

			results['LIGER_validated_signature_collapse_guide'] = expand(resultoutput + 'LIGER/signature-analysis/signature-collapse-guide.tsv')

			analysis_conditions = conditions + ['validated']
			manipulated_conditions = ['validated']

			# if (config['signatures']['signature_collapse_guide'] is not None) & os.path.isfile(config['signatures']['signature_collapse_guide']):
			results['LIGER_collapsed_signature_loading_matrix'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-collapsed.tsv', subtype = ct_extract)
			results['LIGER_collapsed_gene_loading_matrix'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-collapsed.tsv', subtype = ct_extract)

			analysis_conditions = analysis_conditions + ['collapsed']
			manipulated_conditions = manipulated_conditions + ['collapsed']
			report_conditions = ['collapsed']

			if (config['signatures']['score_each_validation_cohort']):
				results['LIGER_collapsed_scored_validation_signature_loading_matrix'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-collapsed-scored-validation.tsv', subtype = ct_extract)
				results['LIGER_collapsed_scored_validation_gene_loading_matrix'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-collapsed-scored-validation.tsv', subtype = ct_extract)

				analysis_conditions = analysis_conditions + ['collapsed-scored-validation']
				manipulated_conditions = manipulated_conditions + ['collapsed-scored-validation']
				report_conditions = report_conditions + ['collapsed-scored-validation']

			results['LIGER_signature_top_gene_loading_matrix_combined'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-loading.tsv', subtype = ct_extract)
			results['LIGER_signature_top_gene_tsv_combined'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene.tsv', subtype = ct_extract)
			results['LIGER_signature_top_gene_loading_heatmap'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-loading.png', subtype = ct_extract)
			results['LIGER_signature_top_gene_loading_matrix_manipulated'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-loading-{condition}.tsv', subtype = ct_extract, condition = manipulated_conditions)
			results['LIGER_signature_top_gene_tsv_manipulated'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-{condition}.tsv', subtype = ct_extract, condition = manipulated_conditions)
			results['LIGER_signature_top_gene_analysis_df_manipulated'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-analysis-{condition}.tsv', subtype = ct_extract, condition = manipulated_conditions)
			results['LIGER_signature_top_gene_loading_heatmap_manipulated'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-loading-{condition}.png', subtype = ct_extract, condition = manipulated_conditions)
			
			results['LIGER_signature_top_gene_analysis_df_manipulated_all_celltype'] = expand(resultoutput + 'LIGER/signature-analysis/signature-top-gene-analysis-combined-{condition}.tsv', condition = ['collapsed', 'collapsed-scored-validation'])

# 			groupedfiles['LIGER_signature_top_markers'] = expand(groupedfilesoutput + 'LIGER/patient-analysis/signature-analysis/top-markers/{condition}/{subtype}/{subtype}-signature-top-gene-loading-{condition}.tsv', subtype = ct_extract, condition = manipulated_conditions)
# 			groupedfiles['LIGER_signature_top_markers_plot'] = expand(groupedfilesoutput + 'LIGER/patient-analysis/signature-analysis/top-markers/{condition}/{subtype}/{subtype}-signature-top-gene-loading-{condition}.png', subtype = ct_extract, condition = manipulated_conditions)

			results['LIGER_signature_geneuniverse'] = []
			results['LIGER_signature_overrepresentation_GO'] = []
			results['LIGER_signature_overrepresentation_KEGG'] = []
			results['LIGER_signature_GSEA_GO'] = []
			results['LIGER_signature_GSEA_3CA'] = []
			results['LIGER_signature_overrepresentation_GO_plot'] = []
			results['LIGER_signature_overrepresentation_KEGG_plot'] = []
			results['LIGER_signature_GSEA_GO_plot'] = []
			results['LIGER_signature_GSEA_3CA_plot'] = []

			results['LIGER_signature_known_markers_loading_matrix_combined'] = []
			results['LIGER_signature_known_markers_gom'] = []
			results['LIGER_signature_known_markers_loading_heatmap'] = []
			results['LIGER_signature_known_markers_gom_heatmap'] = []

			for ct in ct_extract:
				parameterlist_with_idx = parameterlist.set_index('celltype')
				ct_k = parameterlist_with_idx.loc[ct, 'k']
				ct_k_list = list(map(str, np.arange(1, ct_k+1).tolist()))

				results['LIGER_signature_geneuniverse'].extend(expand(resultoutput + 'LIGER/signature-analysis/' + ct + '/enrichment-analysis/' + ct + '-signature-geneuniverse-{condition}.rds', condition = analysis_conditions))
				
				results['LIGER_signature_overrepresentation_GO'].extend(expand(resultoutput + 'LIGER/signature-analysis/' + ct + '/enrichment-analysis/overrepresentation-analysis/{condition}/' + ct + '-signature-{k}-overrepresentation-GO.rds', k = ct_k_list, condition = analysis_conditions))
				results['LIGER_signature_overrepresentation_KEGG'].extend(expand(resultoutput + 'LIGER/signature-analysis/' + ct + '/enrichment-analysis/overrepresentation-analysis/{condition}/' + ct + '-signature-{k}-overrepresentation-KEGG.rds', k = ct_k_list, condition = analysis_conditions))
				results['LIGER_signature_GSEA_GO'].extend(expand(resultoutput + 'LIGER/signature-analysis/' + ct + '/enrichment-analysis/GSEA/{condition}/' + ct + '-signature-{k}-GSEA-GO.rds', k = ct_k_list, condition = analysis_conditions))
				results['LIGER_signature_GSEA_3CA'].extend(expand(resultoutput + 'LIGER/signature-analysis/' + ct + '/enrichment-analysis/GSEA/{condition}/' + ct + '-signature-{k}-GSEA-3CA.rds', k = ct_k_list, condition = analysis_conditions))

				results['LIGER_signature_overrepresentation_GO_plot'].extend(expand(figureoutput + 'LIGER/signature-analysis/' + ct + '/enrichment-analysis/overrepresentation-analysis/{condition}/' + ct + '-signature-{k}-overrepresentation-GO.png', k = ct_k_list, condition = analysis_conditions))
				results['LIGER_signature_overrepresentation_KEGG_plot'].extend(expand(figureoutput + 'LIGER/signature-analysis/' + ct + '/enrichment-analysis/overrepresentation-analysis/{condition}/' + ct + '-signature-{k}-overrepresentation-KEGG.png', k = ct_k_list, condition = analysis_conditions))
				results['LIGER_signature_GSEA_GO_plot'].extend(expand(figureoutput + 'LIGER/signature-analysis/' + ct + '/enrichment-analysis/GSEA/{condition}/' + ct + '-signature-{k}-GSEA-GO.png', k = ct_k_list, condition = analysis_conditions))
				results['LIGER_signature_GSEA_3CA_plot'].extend(expand(figureoutput + 'LIGER/signature-analysis/' + ct + '/enrichment-analysis/GSEA/{condition}/' + ct + '-signature-{k}-GSEA-3CA.png', k = ct_k_list, condition = analysis_conditions))

				if(os.path.isdir(config['signatures']['known_signature_markers_dir'] + ct)):
					if(len(os.listdir(config['signatures']['known_signature_markers_dir'] + ct)) != 0):
						marker_files = os.listdir(config['signatures']['known_signature_markers_dir'] + ct)
						if '.DS_Store' in marker_files:
  							marker_files.remove('.DS_Store')
						marker_lists = [i.split(config['signatures']['known_signature_markers_file_pattern'], 1)[0] for i in marker_files]

						results['LIGER_signature_known_markers_loading_matrix_combined'].extend(expand(resultoutput + 'LIGER/signature-analysis/' + ct + '/gene-loading-analysis/known-markers-loading/' + ct + '-{marker_list}-gene-loading.tsv', marker_list = marker_lists))
						results['LIGER_signature_known_markers_loading_matrix_combined'].extend(expand(resultoutput + 'LIGER/signature-analysis/' + ct + '/gene-loading-analysis/known-markers-loading/' + ct + '-{marker_list}-gene-loading-{condition}.tsv', marker_list = marker_lists, condition = manipulated_conditions))
						results['LIGER_signature_known_markers_gom'].extend(expand(resultoutput + 'LIGER/signature-analysis/' + ct + '/gene-loading-analysis/known-markers-overlap/' + ct + '-{marker_list}-geneoverlap-matrix-{condition}.rds', marker_list = marker_lists, condition = manipulated_conditions))

						results['LIGER_signature_known_markers_loading_heatmap'].extend(expand(figureoutput + 'LIGER/signature-analysis/' + ct + '/gene-loading-analysis/known-markers-loading/' + ct + '-{marker_list}-gene-loading-{filter}.png', marker_list = marker_lists, filter = ['full', 'cleaned']))
						results['LIGER_signature_known_markers_loading_heatmap'].extend(expand(figureoutput + 'LIGER/signature-analysis/' + ct + '/gene-loading-analysis/known-markers-loading/' + ct + '-{marker_list}-gene-loading-{filter}-{condition}.png', marker_list = marker_lists, filter = ['full', 'cleaned'], condition = manipulated_conditions))
						results['LIGER_signature_known_markers_gom_heatmap'].extend(expand(figureoutput + 'LIGER/signature-analysis/' + ct + '/gene-loading-analysis/known-markers-overlap/' + ct + '-{marker_list}-geneoverlap-{stat}-{condition}.png', marker_list = marker_lists, stat = ['oddsratio', 'jaccardindex'], condition = manipulated_conditions))


			if config['signatures']['analyze_sig_loading']:
				results['LIGER_signature_loading_matrix_long_form'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-long-form-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions)
				results['LIGER_signature_loading_uniqueness'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-uniqueness-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions)
				results['LIGER_signature_loading_top_two'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-top-two-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions)
				results['LIGER_signature_top_two_count'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-top-two-count-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions)
				results['LIGER_signature_top_freq'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-top-frequency-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions)
				results['LIGER_signature_loading_patient_summary'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-patient-summary-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions)
				results['LIGER_signature_activation_freq'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-activation-frequency-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions)
				
				plot_grouping = ['grouped', 'summary']

				results['LIGER_signature_loading_qqplot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-qqnorm-{condition}.png', subtype = ct_extract, condition = analysis_conditions)
				results['LIGER_signature_loading_auc_density_grid_plot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-auc-density-grid-{condition}.png', subtype = ct_extract, condition = analysis_conditions)
				results['LIGER_signature_loading_auc_density_grouped_plot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-auc-density-grouped-{condition}.png', subtype = ct_extract, condition = analysis_conditions)
				results['LIGER_signature_loading_auc_box_plot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-loading-auc-box-{grouping}-{condition}.png', subtype = ct_extract, condition = analysis_conditions, grouping = plot_grouping)
				results['LIGER_signature_top_loading_box_plot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-top-loading-box-{grouping}-{condition}.png', subtype = ct_extract, condition = analysis_conditions, grouping = plot_grouping)
				results['LIGER_signature_top_count_box_plot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-top-count-box-{grouping}-{condition}.png', subtype = ct_extract, condition = analysis_conditions, grouping = plot_grouping)
				results['LIGER_signature_top_frequency_bar_plot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/{subtype}-signature-top-freq-bar-{condition}.png', subtype = ct_extract, condition = analysis_conditions)
				
				if config['signatures']['get_sig_loading_profiles'] & (config['signatures']['signature_profile_flavours'] is not None) & os.path.isfile(config['signatures']['signature_profile_flavours']):
					sig_loading_profile_stats = ['mean', 'median']
					sig_profiles = pd.read_csv(config['signatures']['signature_profile_flavours'], header = 0)
					sig_profiles = list(sig_profiles['profile_flavours'])

					results['LIGER_signature_loading_profiles_top_freq'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{condition}/{subtype}-signature-loading-profiles-top-frequency-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions)
					results['LIGER_signature_loading_profiles_loading'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{condition}/{subtype}-signature-loading-profiles-loading-{stat}-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions, stat = sig_loading_profile_stats)
					results['LIGER_signature_loading_profiles_activation_freq'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{condition}/{subtype}-signature-loading-profiles-activation-frequency-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions)

					results['LIGER_signature_loading_profiles_plot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{condition}/{subtype}-signature-loading-profiles-{profile}-{condition}.png', subtype = ct_extract, condition = analysis_conditions, profile = sig_profiles)
					results['LIGER_signature_loading_profiles_corrplot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{condition}/profile-correlations/{subtype}-signature-loading-profiles-{profile}-correlation-{condition}.png', subtype = ct_extract, condition = analysis_conditions, profile = sig_profiles)
					
					results['LIGER_signature_loading_patterns_dimred'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-patterns/{condition}/{subtype}-signature-loading-patterns-dimred-with-top-two-sig-loadings-{condition}.tsv', subtype = ct_extract, condition = report_conditions)

					results['LIGER_signature_top_two_sigs_umap'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-patterns/{condition}/{subtype}-top-two-sigs-umap-{profile}-{condition}.png', subtype = ct_extract, condition = report_conditions, profile = sig_profiles)
					results['LIGER_signature_loading_patterns_umap'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-patterns/{condition}/{subtype}-signature-loading-patterns-umap-{profile}-{condition}.png', subtype = ct_extract, condition = report_conditions, profile = sig_profiles)
					results['LIGER_signature_loading_patterns_corrplot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-patterns/{condition}/pattern-correlations/{subtype}-signature-loading-patterns-{profile}-correlation-{condition}.png', subtype = ct_extract, condition = report_conditions, profile = sig_profiles)
					results['LIGER_signature_loading_patterns_single_cell_corrplot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-patterns/{condition}/pattern-correlations/{subtype}-signature-loading-patterns-single-cell-{profile}-correlation-{condition}.png', subtype = ct_extract, condition = report_conditions, profile = sig_profiles)
					
					results['LIGER_signature_gene_expression'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/gene-expression-analysis/{condition}/{subtype}-signature-gene-expression-{condition}.tsv', subtype = ct_extract, condition = report_conditions)

				if config['signatures']['compute_sig_loading_quantiles']:
					results['LIGER_signature_loading_quantiles_cohort'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-cohort-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions)
					results['LIGER_signature_loading_quantiles_sample'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-sample-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions)
					results['LIGER_signature_loading_quantiles_signature'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-signature-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions)
					results['LIGER_signature_loading_quantiles_signature_cohort'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-signature-cohort-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions)
					results['LIGER_signature_loading_quantiles_signature_sample'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-signature-sample-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions)

					results['LIGER_signature_loading_quantiles_single_var_bar_plot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-single-var-{condition}.png', subtype = ct_extract, condition = analysis_conditions)
					results['LIGER_signature_loading_quantiles_multi_var_box_plot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-analysis/{condition}/signature-quantiles/{subtype}-signature-loading-quantiles-multi-var-{condition}.png', subtype = ct_extract, condition = analysis_conditions)
				
				if (config['patient_profiles']['build_patient_profiles'] & (config['patient_profiles']['profile_ct_config'] is not None) & os.path.isfile(config['patient_profiles']['profile_ct_config'])):
					ct_for_patient_profile = pd.read_csv(config['patient_profiles']['profile_ct_config'], header = 0)
					compartments = list(ct_for_patient_profile.keys())

					patient_profile_flavors = sig_profiles
					
					results['LIGER_patient_signature_profiles'] = expand(resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}.tsv', condition = analysis_conditions, compartment = compartments, profile = patient_profile_flavors)
					results['LIGER_patient_signature_profiles_longform'] = expand(resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-longform-{profile}-{condition}.tsv', condition = analysis_conditions, compartment = compartments, profile = patient_profile_flavors)
					results['LIGER_patient_signature_profiles_correlation'] = expand(resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-correlation-{profile}-{condition}.tsv', condition = analysis_conditions, compartment = compartments, profile = patient_profile_flavors)

					heatmap_groupings = ['grouped', 'clustered']
					
					results['LIGER_patient_signature_profiles_heatmap'] = expand(figureoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}-heatmap-{grouping}.png', condition = analysis_conditions, compartment = compartments, profile = patient_profile_flavors, grouping = heatmap_groupings)
					results['LIGER_patient_signature_profiles_corrplot'] = expand(figureoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}-corrplot.png', condition = analysis_conditions, compartment = compartments, profile = patient_profile_flavors)
					results['LIGER_patient_signature_profiles_stacked_bar_plot'] = expand(figureoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}-stacked-bar-plot.png', condition = analysis_conditions, compartment = compartments, profile = patient_profile_flavors)
					results['LIGER_patient_signature_profiles_circlized_bar_plot'] = expand(figureoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}-circlized-bar-plot.png', condition = analysis_conditions, compartment = compartments, profile = patient_profile_flavors)

					# if (config['signatures']['score_each_validation_cohort']):
					# 	results['LIGER_signature_loading_correlation_collapsed_vs_rescored_validation'] = expand(resultoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-correlation-comparison.tsv', compartment = compartments, profile = patient_profile_flavors)
					# 	results['LIGER_signature_loading_correlation_sign_collapsed_vs_rescored_validation'] = expand(resultoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-correlation-sign-comparison.tsv', compartment = compartments, profile = patient_profile_flavors)
					# 	results['LIGER_signature_loading_correlation_mean_collapsed_vs_rescored_validation'] = expand(resultoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-correlation-mean-comparison.tsv', compartment = compartments, profile = patient_profile_flavors)

					# 	results['LIGER_signature_loading_correlation_scatter_plot_all_intercell_sig_pairs'] = expand(figureoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-all-intercell-sig-pair-correlation-comparison.png', compartment = compartments, profile = patient_profile_flavors)
					# 	results['LIGER_signature_loading_correlation_scatter_plot_intercell_sig_pair_means'] = expand(figureoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-intercell-sig-pair-correlation-mean-comparison.png', compartment = compartments, profile = patient_profile_flavors)
					# 	results['LIGER_signature_loading_correlation_bar_plot_all_sig_pair_corr_sign_comparison'] = expand(figureoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-all-sig-pair-correlation-sign-comparison.png', compartment = compartments, profile = patient_profile_flavors)
					# 	results['LIGER_signature_loading_correlation_bar_plot_sig_pair_corr_same_sign_freq'] = expand(figureoutput + 'LIGER/patient-analysis/collapsed-discovery-vs-rescored-validation/{profile}/patient-{compartment}-signature-loading-{profile}-sig-pair-correlation-same-sign-frequency.png', compartment = compartments, profile = patient_profile_flavors)

				if (config['patient_profiles']['signature_correlation_comparison']):
					results['LIGER_signature_loading_correlation_collapsed_vs_rescored_validation'] = expand(resultoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/patient-{compartment}-signature-profiles-{profile}-correlation-data-frame.tsv', compartment = compartments, profile = patient_profile_flavors)
					results['LIGER_signature_loading_correlation_collapsed_vs_rescored_validation_list'] = expand(resultoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/patient-{compartment}-signature-{profile}-correlation-comparison-data-frame-list.rds', compartment = compartments, profile = patient_profile_flavors)
					results['LIGER_signature_loading_cooccurrence_collapsed_vs_rescored_validation_agreement_data_frame'] = expand(resultoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/patient-{compartment}-signature-{profile}-cooccurrence-agreement-data-frame.tsv', compartment = compartments, profile = patient_profile_flavors)

					results['LIGER_signature_loading_correlation_collapsed_vs_rescored_validation_overall_scatterplot'] = expand(figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/{compartment}-signature-profiles-{profile}-correlation-comparison-overall.png', compartment = compartments, profile = patient_profile_flavors)
					results['LIGER_signature_loading_correlation_collapsed_vs_rescored_validation_intercell_scatterplot'] = expand(figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/{compartment}-signature-profiles-{profile}-correlation-comparison-inter-celltype.png', compartment = compartments, profile = patient_profile_flavors)
					results['LIGER_signature_loading_correlation_collapsed_vs_rescored_validation_intracell_scatterplot'] = expand(figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/{compartment}-signature-profiles-{profile}-correlation-comparison-intra-celltype.png', compartment = compartments, profile = patient_profile_flavors)
					results['LIGER_signature_loading_cooccurrence_collapsed_vs_rescored_validation_agreement_barplot'] = expand(figureoutput + 'LIGER/patient-analysis/signature-correlation-comparison/{profile}/{compartment}-signature-{profile}-cooccurrence-agreement-barplot.png', compartment = compartments, profile = patient_profile_flavors)
				
				if (config['figure2']['make_figure2']):
					results['LIGER_signature_number_df'] = expand(resultoutput + 'LIGER/signature-analysis/figure-2/{compartment}/signature-number-df-{condition}.rds', compartment = compartments, condition = manipulated_conditions)
					results['LIGER_signature_loading_sample_median_mtx'] = expand(resultoutput + 'LIGER/signature-analysis/figure-2/{compartment}/{condition}/signature-loading-median-matrix-{condition}.tsv', compartment = compartments, condition = manipulated_conditions)
					results['LIGER_signature_loading_variance_df'] = expand(resultoutput + 'LIGER/signature-analysis/figure-2/{compartment}/{condition}/signature-loading-variance-{condition}.tsv', compartment = compartments, condition = manipulated_conditions)
					results['LIGER_signature_top_gene_df_list'] = expand(resultoutput + 'LIGER/signature-analysis/figure-2/{compartment}/{condition}/selected-signature-top-gene-loading-mtrices-{condition}.rds', compartment = compartments, condition = manipulated_conditions)

					results['figure_2_signature_number'] = expand(figureoutput + 'LIGER/signature-analysis/figure-2/{compartment}/figure-2-a-{compartment}.png', compartment = compartments)
					results['figure_2_signature_heatmap'] = expand(figureoutput + 'LIGER/signature-analysis/figure-2/{compartment}/figure-2-b-{compartment}.png', compartment = compartments)
					results['figure_2_signature_markers'] = expand(figureoutput + 'LIGER/signature-analysis/figure-2/{compartment}/figure-2-c-{compartment}.png', compartment = compartments)
					results['figure_2_plot'] = expand(figureoutput + 'LIGER/signature-analysis/figure-2/{compartment}/figure-2-{compartment}.png', compartment = compartments)

				if (config['figure1']['make_figure1']):
					results['figure_1_new_plot'] = expand(figureoutput + 'LIGER/signature-analysis/figure-1-new/figure-1.png')

				if (config['figure3']['make_figure3']):
					results['LIGER_signature_loading_variance_df_combined'] = expand(resultoutput + 'LIGER/signature-analysis/figure-3/{compartment}/signature-loading-variance-dis-val.tsv', compartment = compartments)

					results['figure_3_signature_loading_variance_combined'] = expand(figureoutput + 'LIGER/signature-analysis/figure-3/{compartment}/figure-3-a-{compartment}.png', compartment = compartments)
					results['figure_3_signature_loading_variance_ct_facet'] = expand(figureoutput + 'LIGER/signature-analysis/figure-3/{compartment}/figure-3-b-{compartment}.png', compartment = compartments)
					results['figure_3_plot'] = expand(figureoutput + 'LIGER/signature-analysis/figure-3/{compartment}/figure-3-{compartment}.png', compartment = compartments)

				if (config['figure4']['make_figure4']):
					results['LIGER_signature_loading_correlation_df_full_and_intra'] = expand(resultoutput + 'LIGER/patient-analysis/figure-4/{profile}/patient-{compartment}-signature-profiles-{profile}-correlation-data-frame-full-and-intra.tsv', compartment = compartments, profile = patient_profile_flavors)
					results['LIGER_signature_loading_correlation_df_inter'] = expand(resultoutput + 'LIGER/patient-analysis/figure-4/{profile}/patient-{compartment}-signature-profiles-{profile}-correlation-data-frame-inter.tsv', compartment = compartments, profile = patient_profile_flavors)
					results['LIGER_signature_loading_cooccurrence_agreement_examples'] = expand(resultoutput + 'LIGER/patient-analysis/figure-4/{profile}/patient-{compartment}-signature-{profile}-correlation-comparison-data-frame-list-examples.rds', compartment = compartments, profile = patient_profile_flavors)

					results['figure_4_plot'] = expand(figureoutput + 'LIGER/patient-analysis/figure-4/{compartment}/loading-mean/figure-4-loading-mean-{compartment}.png', compartment = compartments)
				
				if (config['figure4']['make_figure4']):
					results['figure_2_new_plot'] = expand(figureoutput + 'LIGER/patient-analysis/figure-2-new/{compartment}/loading-mean/figure-2-loading-mean-{compartment}.png', compartment = compartments)
				
				if(config['stan']['run_stan'] & (config['stan']['stan_ct_config'] is not None) & os.path.isfile(config['stan']['stan_ct_config']) & os.path.isfile(config['stan']['stan_param_config'])):
					stan_ct_config = pd.read_csv(config['stan']['stan_ct_config'], header = 0)
					scopes = list(stan_ct_config.keys())

					stan_parameters = pd.read_csv(config['stan']['stan_param_config'], header = 0)
					num_niches_list = list(stan_parameters['num_niches'])
					nIter_list = list(stan_parameters['nIter'])

					stan_conditions = ['collapsed'] + ['collapsed-scored-validation']

					results['stan_data'] = expand(resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-data-{scope}-env-{condition}-num-niches-{numniche}.rds', scope = scopes, condition = stan_conditions, numniche = num_niches_list)
					results['stan_df_sig_mean'] = expand(resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-df-sig-mean-{scope}-env-{condition}-num-niches-{numniche}.rds', scope = scopes, condition = stan_conditions, numniche = num_niches_list)
					results['stan_samples_encodings'] = expand(resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-samples-encodings-{scope}-env-{condition}-num-niches-{numniche}.tsv', scope = scopes, condition = stan_conditions, numniche = num_niches_list)
					results['stan_sigs_encodings'] = expand(resultoutput + 'LIGER/patient-analysis/stan/data/{condition}/stan-sigs-encodings-{scope}-env-{condition}-num-niches-{numniche}.tsv', scope = scopes, condition = stan_conditions, numniche = num_niches_list)
					results['stan_model'] = expand(resultoutput + 'LIGER/patient-analysis/stan/model/model-{scope}-env-{condition}-num-niches-{numniche}-' + config['stan']['stan_model_basename'] + '.stan', scope = scopes, condition = stan_conditions, numniche = num_niches_list)

					results['fit_optim'] = expand(resultoutput + 'LIGER/patient-analysis/stan/results/model-fit/{condition}/stan-model-fit-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds', scope = scopes, condition = stan_conditions, numniche = num_niches_list, niter = nIter_list)

					results['patient_specific_modelled_mu'] = expand(resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/patient-specific-modelled-mu-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds', scope = scopes, condition = stan_conditions, numniche = num_niches_list, niter = nIter_list)
					results['microenvironment_niche_factors'] = expand(resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/microenvironment-niche-factors-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds', scope = scopes, condition = stan_conditions, numniche = num_niches_list, niter = nIter_list)
					results['niche_factor_loadings'] = expand(resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/niche-factor-loadings-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds', scope = scopes, condition = stan_conditions, numniche = num_niches_list, niter = nIter_list)
					results['intrinsic_covariance_matrices'] = expand(resultoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/intrinsic-covariance-matrices-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.rds', scope = scopes, condition = stan_conditions, numniche = num_niches_list, niter = nIter_list)

					results['patient_specific_modelled_mu_plot'] = expand(figureoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/patient-specific-modelled-mu-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.png', scope = scopes, condition = stan_conditions, numniche = num_niches_list, niter = nIter_list)
					results['microenvironment_niche_factors_plot'] = expand(figureoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/microenvironment-niche-factors-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.png', scope = scopes, condition = stan_conditions, numniche = num_niches_list, niter = nIter_list)
					results['niche_factor_loadings_plot'] = expand(figureoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/niche-factor-loadings-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.png', scope = scopes, condition = stan_conditions, numniche = num_niches_list, niter = nIter_list)
					results['intrinsic_covariance_matrices_plot'] = expand(figureoutput + 'LIGER/patient-analysis/stan/results/parameters/{condition}/intrinsic-covariance-matrices-{scope}-env-{condition}-num-niches-{numniche}-niter-{niter}.png', scope = scopes, condition = stan_conditions, numniche = num_niches_list, niter = nIter_list)
					results['microenvironment_niche_factors_plot_combined'] = expand(figureoutput + 'LIGER/patient-analysis/stan/results/parameter-comparison/microenvironment-niche-factors-{scope}-env-collapsed-and-scored-validation-num-niches-{numniche}-niter-{niter}-sig-loading.png', scope = scopes, numniche = num_niches_list, niter = nIter_list)
					results['microenvironment_niche_factors_correlation_plot'] = expand(figureoutput + 'LIGER/patient-analysis/stan/results/parameter-comparison/microenvironment-niche-factors-{scope}-env-collapsed-vs-scored-validation-num-niches-{numniche}-niter-{niter}-correlation.png', scope = scopes, numniche = num_niches_list, niter = nIter_list)
					results['intrinsic_covariance_matrices_correlation_plot'] = expand(figureoutput + 'LIGER/patient-analysis/stan/results/parameter-comparison/intrinsic-covariance-matrices-{scope}-env-collapsed-vs-scored-validation-num-niches-{numniche}-niter-{niter}-correlation.png', scope = scopes, numniche = num_niches_list, niter = nIter_list)
					results['intrinsic_covariance_matrices_correlation_bar_plot'] = expand(figureoutput + 'LIGER/patient-analysis/stan/results/parameter-comparison/intrinsic-covariance-matrices-{scope}-env-collapsed-vs-scored-validation-num-niches-{numniche}-niter-{niter}-correlation-bar.png', scope = scopes, numniche = num_niches_list, niter = nIter_list)

					if(config['figure5']['make_figure5']):
						results['microenvironment_niche_factors_loading_combined'] = expand(resultoutput + 'LIGER/patient-analysis/figure-5/microenvironment-niche-factors-{scope}-env-combined.tsv', scope = scopes)
						results['intrinsic_covariance_matrices_correlation_data'] = expand(resultoutput + 'LIGER/patient-analysis/figure-5/intrinsic-covariance-matrices-{scope}-env-corr.tsv', scope = scopes)
						results['intrinsic_covariance_matrices_list_to_plot'] = expand(resultoutput + 'LIGER/patient-analysis/figure-5/intrinsic-covariance-matrices-{scope}-env.rds', scope = scopes)

						results['figure_5_plot'] = expand(figureoutput + 'LIGER/patient-analysis/figure-5/{scope}/figure-5-{scope}.png', scope = scopes)

					if(config['figure5']['make_figure5']):
						results['figure_3_new_plot'] = expand(figureoutput + 'LIGER/patient-analysis/figure-3-new/{scope}/figure-3-{scope}.png', scope = scopes)

					if(config['survival']['run_survival_analysis']):
						results['survival_analysis_TCGA_PDAC_clinical_data'] = expand(resultoutput + 'LIGER/survival-analysis/clinical-data/paad-clin-data.tsv')
						results['survival_analysis_COMPASS_clinical_data'] = expand(resultoutput + 'LIGER/survival-analysis/clinical-data/compass-clin-data.tsv')
						results['survival_analysis_PanCuRx_clinical_data'] = expand(resultoutput + 'LIGER/survival-analysis/clinical-data/pancurx-clin-data.tsv')
						results['survival_analysis_Toronto_COMPASS_clin_data_comparison'] = expand(resultoutput + 'LIGER/survival-analysis/clinical-data/toronto-compass-clin-data-comparison.tsv')
						
						results['survival_analysis_TCGA_PDAC_expression_data'] = expand(resultoutput + 'LIGER/survival-analysis/bulk-rna-data/paad-mrna-data.rds')
						results['survival_analysis_COMPASS_expression_data'] = expand(resultoutput + 'LIGER/survival-analysis/bulk-rna-data/compass-mrna-data.rds')
						results['survival_analysis_PanCuRx_expression_data'] = expand(resultoutput + 'LIGER/survival-analysis/bulk-rna-data/pancurx-mrna-data.rds')
						results['survival_analysis_COMPASS_Pa_expression_data'] = expand(resultoutput + 'LIGER/survival-analysis/bulk-rna-data/compass-pa-mrna-data.rds')
						results['survival_analysis_COMPASS_Lv_expression_data'] = expand(resultoutput + 'LIGER/survival-analysis/bulk-rna-data/compass-lv-mrna-data.rds')

						results['survival_analysis_signature_TCGA_PDAC_CM_curve_plot_list'] = expand(resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/cm-curve-plot-list-paad.rds', scope = scopes)
						results['survival_analysis_signature_PanCuRx_CM_curve_plot_list'] = expand(resultoutput + 'LIGER/survival-analysis/{scope}/signature-survival-analysis/cm-curve-plot-list-pancurx.rds', scope = scopes)
						results['survival_analysis_niche_CM_curve_plot_list'] = expand(resultoutput + 'LIGER/survival-analysis/{scope}/niche-survival-analysis/cm-curve-plot-list-paad-niche.rds', scope = scopes)
						results['survival_analysis_niche_PanCuRx_CM_curve_plot_list'] = expand(resultoutput + 'LIGER/survival-analysis/{scope}/niche-survival-analysis/cm-curve-plot-list-pancurx-niche.rds', scope = scopes)

						results['figure_6_plot'] = expand(figureoutput + 'LIGER/survival-analysis/figure-6/figure-6.png')
						results['figure_4_new_plot'] = expand(figureoutput + 'LIGER/survival-analysis/figure-4-new/figure-4.png')

						results['figure_niche_loadings_in_bulkrna'] = expand(figureoutput + 'LIGER/survival-analysis/full/niche-survival-analysis/niche-loadings-pancurx.png')
					

			

### [ IMPORT INDIVIDUAL SNAKEMAKE FILES ] #####
# import snakemake sub-workflows
include: 'pipeline/process-data.smk'
include: 'pipeline/cell-type-assignment.smk'
include: 'pipeline/cell-selection.smk'
include: 'pipeline/data-organization.smk'
include: 'pipeline/signature-extraction.smk'
include: 'pipeline/signature-validation.smk'
include: 'pipeline/signature-analysis.smk'
include: 'pipeline/patient-analysis.smk'
include: 'pipeline/survival-analysis.smk'
# include: 'pipeline/output-files-grouping.smk'

# Save report as
report: "report/workflow.rst"

rule all:
    input:
        results.values()#,
	#groupedfiles.values(),
	# report
