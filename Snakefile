### [ SET UP ] #####
import pandas as pd
import numpy as np 
import yaml
from itertools import chain
import os

configfile: 'config/config.yml'

# Functions 
def get_mem_mb(attempt, base_mem):
    return base_mem + (attempt * 2000)

# Input files can be passed as references to the configuration file, however, this only works when the config file is not empty
def get_config_files(param):
    config[param]


cohorts = config['cohorts']

# This determines the output folder - the version can be set in the config file
output = 'output/' + config['version'] + '/'
figureoutput = output + 'figures/'
dataoutput = output + 'data/'
resultoutput = output + 'results/'

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

def get_liger_param(parameterlist, celltype, param):
	ct_parameterlist = parameterlist[parameterlist['celltype'] == celltype]
	return list(ct_parameterlist[param])[0]

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
# outputs
results = {
	'combined_doublet_detection': expand(resultoutput + 'doublet-detection/DoubletFinder-combined-{cohort}.tsv', cohort = cohorts),

    	'scRNA_filtered_sce': expand(dataoutput + 'process-sce/scRNASeq-filtered-sce-{cohort}.rds', cohort = cohorts),
	'cells_per_sample': expand(resultoutput + 'process-sce/number-of-cells-{cohort}.tsv', cohort = cohorts),
    	'cells_per_sample_plot': expand(figureoutput + 'process-sce/number-of-cells-{cohort}.png', cohort = cohorts),

	'scRNA_Azimuth_assigned_seu': expand(dataoutput + 'Azimuth-annotation/scRNASeq-Azimuth-assigned-seu-{cohort}.rds', cohort = cohorts),
    	'scRNA_Azimuth_assigned_sce': expand(dataoutput + 'Azimuth-annotation/scRNASeq-Azimuth-assigned-sce-{cohort}.rds', cohort = cohorts),
	'cell_type_assignment_scores': expand(resultoutput + 'Azimuth-annotation/cell-type-annotation-scores-{cohort}.tsv', cohort = cohorts),

	'cell_type_and_sample_UMAP_plot': expand(figureoutput + 'Azimuth-annotation/annotated-dimred/UMAP-cell-type-and-sample-{cohort}.png', cohort = cohorts),
	'cell_type_prediction_score_UMAP_plot': expand(figureoutput + 'Azimuth-annotation/celltype-prediction-score/UMAP-celltype-prediction-score-{cohort}.png', cohort = cohorts),
	'cell_type_prediction_score_heatmap': expand(figureoutput + 'Azimuth-annotation/celltype-prediction-score/heatmap-celltype-prediction-score-{cohort}.png', cohort = cohorts),
}

if (config['visualize_these_cohorts_together'] is not None) & os.path.isfile(config['visualize_these_cohorts_together']): 
	# Visualize cell type labels predicted by Azimuth for all cohorts in one plot
	cohorts_together = pd.read_csv(config['visualize_these_cohorts_together'], header = 0)
	cohorts_together = list(cohorts_together['cohort'])

	colour_field = ['cohort', 'celltype']

	results['combined_UMAP_plots'] = expand(figureoutput + 'Azimuth-annotation/annotated-dimred/cohort-combined-plots/UMAP-{field}-all-cohorts.png', field = colour_field) 

if unique_markers:
	results['cell_type_marker_plot'] = expand(figureoutput + 'Azimuth-annotation/celltype-markers/{cohort}/{marker}.png', cohort = cohorts, marker = unique_markers)

if (config['celltypes_csv'] is not None) & os.path.isfile(config['celltypes_csv']):
	# Create a list of all cell types of interest
	celltypes = pd.read_csv(config['celltypes_csv'], header = 0) 
	celltypes = list(celltypes['celltype'])
	
	celltype_hierarchy = {k: [] for k in celltypes}

	results['scRNA_subsetted_sce'] = expand(dataoutput + 'subset-sce/{celltype}/scRNASeq-{celltype}-sce-{cohort}.rds', cohort = cohorts, celltype = celltypes)
	results['subsetted_cells_per_sample_plot'] = expand(figureoutput + 'subset-sce/{celltype}/number-of-{celltype}-cells-{cohort}.png', cohort = cohorts, celltype = celltypes)
	results['subsetted_cells_per_sample'] = expand(resultoutput + 'subset-sce/{celltype}/number-of-{celltype}-cells-{cohort}.tsv', cohort = cohorts, celltype = celltypes)

if (config['infercnv']['gene_order_file'] is not None) & os.path.isfile(config['infercnv']['gene_order_file']):
	gene_order_file_path = config['infercnv']['gene_order_file']
	obs_cohorts = pd.read_csv(config['infercnv']['obs_cohorts'], header = 0)
	obs_cohorts = list(obs_cohorts['cohort'])

	results['inferCNV_objects'] = expand(resultoutput + 'inferCNV/' + config['infercnv']['celltype'] + '/{cohort}/ref-' + config['infercnv']['ref_cohort'] + '/run.final.infercnv_obj', cohort = obs_cohorts)

if (config['singler']['celltypes_focused'] is not None) & os.path.isfile(config['singler']['celltypes_focused']):
	# Create a list of cell types for subtype annotation
	singler = pd.read_csv(config['singler']['celltypes_focused'], header = 0)
	singler_info = dict(zip(singler['celltype'], singler['ref_label_field']))
	ct_zoomin = singler_info.keys()

	results['scRNA_SingleR_annotated_sce'] = expand(dataoutput + 'SingleR-annotation/{celltype}/sceRNASeq-SingleR-annotated-{celltype}-sce-{cohort}.rds', cohort = cohorts, celltype = ct_zoomin)
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

		if ct == "immune":
			results['immune_cell_type_marker_plot'] = expand(figureoutput + 'SingleR-annotation/{celltype}/celltype-markers/{cohort}/{marker}.png', cohort = cohorts, celltype = ct, marker = unique_immune_markers),

if os.path.isfile(config['celltype_annot_resource_dir'] + 'cell_types_to_combine.csv'):
	ct_zoomout_info = pd.read_csv(config['celltype_annot_resource_dir'] + 'cell_types_to_combine.csv', header = 0)
	ct_zoomout = list(ct_zoomout_info.keys())
	results['scRNA_combined_sce'] = []

	# for ct in ct_zoomout:
	# 	results['scRNA_combined_sce'].extend(expand(get_combined_sce_dir(ct, ct_zoomout_info, dataoutput, celltype_hierarchy) + ct + '/scRNASeq-' + ct + '-sce-{cohort}.rds', cohort = cohorts))

if (config['liger']['celltypes_for_signature'] is not None) & os.path.isfile(config['liger']['celltypes_for_signature']):
	# Create a list of cell types for extracting signatures
	ct_signature = pd.read_csv(config['liger']['celltypes_for_signature'], header = 0)
	ct_signature = list(ct_signature['celltype'])
	
	if (config['liger']['discovery_cohort_list'] is not None) & os.path.isfile(config['liger']['discovery_cohort_list']):
		discovery_cohorts = pd.read_csv(config['liger']['discovery_cohort_list'], header = 0)
		discovery_cohorts = list(discovery_cohorts['cohort'])
		validation_cohorts = pd.read_csv(config['liger']['validation_cohort_list'], header = 0)
		validation_cohorts = list(validation_cohorts['cohort'])

		results['scRNA_discovery_sce_combined'] = expand(dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNAseq-{subtype}-sce-discovery.rds', subtype = ct_signature)
		results['scRNA_validation_sce_combined'] = expand(dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNAseq-{subtype}-sce-validation.rds', subtype = ct_signature)
		results['scRNA_discovery_sce_list'] = expand(dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNAseq-{subtype}-scelist-discovery.rds', subtype = ct_signature)
		results['scRNA_validation_sce_list'] = expand(dataoutput + 'cohort-discovery-validation-grouping/{subtype}/scRNAseq-{subtype}-scelist-validation.rds', subtype = ct_signature)
		
		results['discovery_cells_per_sample'] = expand(resultoutput + 'cohort-discovery-validation-grouping/{subtype}/number-of-{subtype}-cells-discovery.tsv', subtype = ct_signature)
		results['validation_cells_per_sample'] = expand(resultoutput + 'cohort-discovery-validation-grouping/{subtype}/number-of-{subtype}-cells-validation.tsv', subtype = ct_signature)

		if (config['liger']['seedlist'] is not None) & os.path.isfile(config['liger']['seedlist']):
			seeds = pd.read_csv(config['liger']['seedlist'], header = 0)
			seeds = list(seeds['seed'])
			conditions = ['discovery', 'validation']

			if (config['liger']['parameter_range'] is not None) & os.path.isfile(config['liger']['parameter_range']):
				parameter_range = pd.read_csv(config['liger']['parameter_range'], header = 0)
				lambda_range = list(parameter_range['lambda'])
				k_range = list(parameter_range['k'])
				lambdalist = list(range(lambda_range[0], lambda_range[1]))
				klist = list(range(k_range[0], k_range[1]))
			else:
				lambdalist = np.arange(1, 5, 1).tolist()
				klist = np.arange(5, 20, 5).tolist()

			results['LIGER_parameter_sweep_metrics'] = expand(resultoutput + 'LIGER/parameter-sweep/metrics/{subtype}/{condition}/seed-{seed}/{subtype}-{condition}-seed-{seed}-lambda-{Lambda}-k-{k}.tsv', subtype = ct_signature, condition = conditions, seed = seeds, Lambda = lambdalist, k = klist)
			results['LIGER_parameter_sweep_factor_contributions'] = expand(resultoutput + 'LIGER/parameter-sweep/factor-contributions/{subtype}/{condition}/seed-{seed}/k-{k}/{subtype}-{condition}-seed-{seed}-lambda-{Lambda}-k-{k}.tsv', subtype = ct_signature, condition = conditions, seed = seeds, Lambda = lambdalist, k = klist)

			results['LIGER_parameter_sweep_metrics_combined'] = expand(resultoutput + 'LIGER/parameter-sweep/metrics/{subtype}/{condition}/{subtype}-{condition}-parameter-sweep-metrics.tsv', subtype = ct_signature, condition = conditions)
			results['LIGER_parameter_sweep_metrics_plot_k'] = expand(figureoutput + 'LIGER/parameter-sweep/metrics/{subtype}/{condition}/{subtype}-{condition}-k-sweep.png', subtype = ct_signature, condition = conditions)
			results['LIGER_parameter_sweep_metrics_plot_lambda'] = expand(figureoutput + 'LIGER/parameter-sweep/metrics/{subtype}/{condition}/{subtype}-{condition}-lambda-sweep.png', subtype = ct_signature, condition = conditions)

		if (config['liger']['parameterlist'] is not None) & os.path.isfile(config['liger']['parameterlist']):
			parameterlist = pd.read_csv(config['liger']['parameterlist'], header = 0)
			ct_extract = list(parameterlist['celltype'])
			
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

			results['LIGER_validated_signature_loading_matrix'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-validated.tsv', subtype = ct_extract)
			results['LIGER_validated_gene_loading_matrix'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-validated.tsv', subtype = ct_extract)

			analysis_conditions = conditions + ['validated']
			manipulated_conditions = ['validated']

			if (config['signatures']['signature_collapse_guide'] is not None) & os.path.isfile(config['signatures']['signature_collapse_guide']):
				results['LIGER_collapsed_signature_loading_matrix'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-signature-loading-collapsed.tsv', subtype = ct_extract)
				results['LIGER_collapsed_gene_loading_matrix'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/loading-matrices/{subtype}-gene-loading-collapsed.tsv', subtype = ct_extract)

				analysis_conditions = analysis_conditions + ['collapsed']
				manipulated_conditions = manipulated_conditions + ['collapsed']

			results['LIGER_signature_top_gene_loading_matrix_combined'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-loading.tsv', subtype = ct_extract)
			results['LIGER_signature_top_gene_loading_heatmap'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-loading.png', subtype = ct_extract)
			results['LIGER_signature_top_gene_loading_matrix_combined_manipulated'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-loading-{condition}.tsv', subtype = ct_extract, condition = manipulated_conditions)
			results['LIGER_signature_top_gene_loading_heatmap_manipulated'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-loading-{condition}.png', subtype = ct_extract, condition = manipulated_conditions)

			results['LIGER_signature_geneuniverse'] = []
			results['LIGER_signature_overrepresentation_GO'] = []
			results['LIGER_signature_overrepresentation_KEGG'] = []
			results['LIGER_signature_GSEA_GO'] = []
			results['LIGER_signature_overrepresentation_GO_plot'] = []
			results['LIGER_signature_overrepresentation_KEGG_plot'] = []
			results['LIGER_signature_GSEA_GO_plot'] = []

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

				results['LIGER_signature_overrepresentation_GO_plot'].extend(expand(figureoutput + 'LIGER/signature-analysis/' + ct + '/enrichment-analysis/overrepresentation-analysis/{condition}/' + ct + '-signature-{k}-overrepresentation-GO.png', k = ct_k_list, condition = analysis_conditions))
				results['LIGER_signature_overrepresentation_KEGG_plot'].extend(expand(figureoutput + 'LIGER/signature-analysis/' + ct + '/enrichment-analysis/overrepresentation-analysis/{condition}/' + ct + '-signature-{k}-overrepresentation-KEGG.png', k = ct_k_list, condition = analysis_conditions))
				results['LIGER_signature_GSEA_GO_plot'].extend(expand(figureoutput + 'LIGER/signature-analysis/' + ct + '/enrichment-analysis/GSEA/{condition}/' + ct + '-signature-{k}-GSEA-GO.png', k = ct_k_list, condition = analysis_conditions))

				if(os.path.isdir(config['signatures']['known_signature_markers_dir'] + ct)):
					if(len(os.listdir(config['signatures']['known_signature_markers_dir'] + ct)) != 0):
						marker_files = os.listdir(config['signatures']['known_signature_markers_dir'] + ct)
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
				
				if config['signatures']['get_sig_loading_profiles']:
					sig_loading_profile_stats = ['mean', 'median']

					results['LIGER_signature_loading_profiles_top_freq'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{condition}/{subtype}-signature-loading-profiles-top-frequency-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions)
					results['LIGER_signature_loading_profiles_loading'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{condition}/{subtype}-signature-loading-profiles-loading-{stat}-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions, stat = sig_loading_profile_stats)
					results['LIGER_signature_loading_profiles_activation_freq'] = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{condition}/{subtype}-signature-loading-profiles-activation-frequency-{condition}.tsv', subtype = ct_extract, condition = analysis_conditions)

					sig_profiles = ['top-frequency', 'loading-mean', 'loading-median', 'activation-frequency']

					results['LIGER_signature_loafing_profiles_plot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{condition}/{subtype}-signature-loading-profiles-{profile}-{condition}.png', subtype = ct_extract, condition = analysis_conditions, profile = sig_profiles)
					results['LIGER_signature_loafing_profiles_corrplot'] = expand(figureoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{condition}/profile-correlations/{subtype}-signature-loading-profiles-{profile}-correlation-{condition}.png', subtype = ct_extract, condition = analysis_conditions, profile = sig_profiles)

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

			

### [ IMPORT INDIVIDUAL SNAKEMAKE FILES ] #####
# import snakemake sub-workflows
include: 'pipeline/process-data.smk'
include: 'pipeline/cell-type-assignment.smk'
include: 'pipeline/cell-selection.smk'
include: 'pipeline/data-organization.smk'
include: 'pipeline/signature-extraction.smk'
include: 'pipeline/signature-analysis.smk'
include: 'pipeline/patient-analysis.smk'

# Save report as
report: "report/workflow.rst"

rule all:
    input:
        results.values()
