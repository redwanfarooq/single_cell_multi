##########################################################################################
# Snakemake rule for integrating multisample multimodal data
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
# Requires outputs from resources/rules/normalise.smk
##########################################################################################

scripts_dir = config.get("scripts_dir", "resources/scripts")

# Define rule
rule integrate:
	input: normalise,
	output: os.path.join(config["output_dir"], "integrate", f"{{integrate}}.{config.get('format', 'qs')}"), 
	log: os.path.abspath("logs/integrate/{integrate}.log")
	threads: threads: min(round(len(samples) / 5), 4)
	resources:
		mem = lambda wildcards, threads: f"{threads * 100}GiB"
	params:
		script_path = scripts_dir if os.path.isabs(scripts_dir) else os.path.join(workflow.basedir, scripts_dir),
		optional_flags = get_optional_flags(rna_assay = config.get("rna-assay", None), atac_assay = config.get("atac-assay", None), adt_assay = config.get("adt-assay", None), batch=";".join(config.get("batch", list())), normalisation_method=config.get("normalisation-method", None), integration_method=parse_integration_method(config.get("integration-method", None)), integration_args=config.get("integration-args", None)),
	conda: "single_cell_multi"
	# envmodules: "R-cbrg"
	message: "Integrating data"
	shell:
		"""
		cd {params.script_path} && \
		./integrate.R \
			--input {input} \
			{params.optional_flags} \
			--output {output} \
			--threads {threads} \
			--log {log}
		"""

	
# Set rule targets
integrate = expand(os.path.join(config["output_dir"], "integrate", f"batch:{{batch}}_normalisation:{{normalisation_method}}_integration:{{integration_method}}.{config.get('format', 'qs')}"), batch="+".join(config.get("batch", ["orig.ident"])), normalisation_method=config.get("normalisation-method", "LogNormalize"), integration_method=config.get("integration-method", "rpca"))