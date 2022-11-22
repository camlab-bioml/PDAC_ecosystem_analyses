rule construct_liger_signature_patient_profiles:
	input:
		#celltype_sig_profiles = expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/{{condition}}/{subtype}-signature-loading-profiles-{{profile}}-{{condition}}.tsv', subtype = ),
		celltype_sig_profiles = lambda wildcards: expand(resultoutput + 'LIGER/signature-analysis/{subtype}/signature-loading-profiles/' + wildcards.condition + '/{subtype}-signature-loading-profiles-' + wildcards.profile + '-' + wildcards.condition + '.tsv', subtype = list(ct_for_patient_profile[wildcards.compartment])),
	params:
	output:
		patient_profiles = resultoutput + 'LIGER/patient-analysis/patient-signature-profiles/{condition}/{profile}/patient-{compartment}-signature-profiles-{profile}-{condition}.tsv',
	resources:	
		mem_mb = 1000
	threads: 1
	script: 
		'patient-analysis/construct-liger-signature-patient-profiles.R'
