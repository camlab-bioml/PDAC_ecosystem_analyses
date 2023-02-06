rule construct_liger_signature_patient_profiles:
	input:
		celltype_sig_profiles = lambda wildcards: expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/' + wildcards.condition + '/{subtype}-signature-loading-profiles-' + wildcards.profile + '-' + wildcards.condition + '.tsv', subtype = list(ct_for_patient_profile[wildcards.compartment])),
	params:
	output:
		patient_profiles = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}.tsv',
		patient_profiles_long = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-longform-{profile}-{condition}.tsv',
		patient_profiles_corr = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-correlation-{profile}-{condition}.tsv',
	resources:	
		mem_mb = 1000
	threads: 1
	script: 
		'patient-analysis/construct-liger-signature-patient-profiles.R'

rule visualize_liger_signature_patient_profiles:
	input:
		patient_profiles = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}.tsv',
		patient_profiles_long = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-longform-{profile}-{condition}.tsv',
		patient_profiles_corr = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-correlation-{profile}-{condition}.tsv',
	params:
	output:
		patient_profiles_heatmap_grouped = figureoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}-heatmap-grouped.png',
		patient_profiles_heatmap_clustered = figureoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}-heatmap-clustered.png',
		patient_profiles_corrplot = figureoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}-corrplot.png',
		patient_profiles_stacked_bar_plot = figureoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}-stacked-bar-plot.png',
		patient_profiles_circlized_bar_plot = figureoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}-circlized-bar-plot.png',
	resources:	
		mem_mb = 1000
	threads: 1
	script: 
		'patient-analysis/visualize-liger-signature-patient-profiles.R'
