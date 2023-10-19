rule group_output_files_liger_signature_top_markers:
	input:
		sig_top_gene_loading_mtx = resultoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-loading-{condition}.tsv',
		sig_top_gene_loading_heatmap = figureoutput + 'LIGER/signature-analysis/{subtype}/gene-loading-analysis/{subtype}-signature-top-gene-loading-{condition}.png',
	params:
	output:
		sig_top_gene_loading_mtx = groupedfilesoutput + 'LIGER/patient-analysis/signature-analysis/top-markers/{condition}/{subtype}/{subtype}-signature-top-gene-loading-{condition}.tsv',
		sig_top_gene_loading_heatmap = groupedfilesoutput + 'LIGER/patient-analysis/signature-analysis/top-markers/{condition}/{subtype}/{subtype}-signature-top-gene-loading-{condition}.png',

	resources:	
		mem_mb = 1000
	threads: 1
	shell: 
		'cp {input.sig_top_gene_loading_mtx} {output.sig_top_gene_loading_mtx} && cp {input.sig_top_gene_loading_heatmap} {output.sig_top_gene_loading_heatmap}'