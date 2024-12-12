##########################################################################################
# Snakemake rule for integrating multisample multimodal data
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
# Requires outputs from resources/rules/normalise.smk
##########################################################################################

scripts_dir = config.get("scripts_dir", "resources/scripts")

# Define rule
rule integrate:
	input: os.path.join(config["output_dir"], "{prefix}", f"normalised.{config.get('format', 'qs')}"),
	output: os.path.join(config["output_dir"], "{prefix}", "batch:{batch}_normalisation:{normalisation_method}_integration:{integration_method}", f"integrated.{config.get('format', 'qs')}"), 
	log: os.path.abspath("logs/{prefix}/batch:{batch}_normalisation:{normalisation_method}_integration:{integration_method}/integrate.log")
	threads: lambda wildcards: min(floor(len(samples) / 5) + 1, 4) if wildcards.integration_method != "scvi" else 1
	resources:
		partition = lambda wildcards: "gpu" if wildcards.integration_method == "scvi" else "short",
		gpus = lambda wildcards: 1 if wildcards.integration_method == "scvi" else 0,
		mem = lambda wildcards, threads: f"{threads * 100}GiB" if wildcards.integration_method != "scvi" else f"{min(len(samples) * 5, 96)}GiB",
	params:
		script_path = scripts_dir if os.path.isabs(scripts_dir) else os.path.join(workflow.basedir, scripts_dir),
		wildcard_flags = lambda wildcards: get_optional_flags(batch=";".join(wildcards.batch.split("+")), normalisation_method=wildcards.normalisation_method, integration_method=parse_integration_method(wildcards.integration_method) if wildcards.integration_method != "None" else None, integration_args=config.get("integration-args", {wildcards.integration_method: None}).get(wildcards.integration_method, None)),
		other_flags = get_optional_flags(rna_assay = config.get("rna-assay", None), atac_assay = config.get("atac-assay", None), adt_assay = config.get("adt-assay", None)),
	conda: "single_cell_multi"
	# envmodules: "R-cbrg"
	message: "Performing dimensionality reduction and integrating data"
	shell:
		"""
		cd {params.script_path} && \
		./integrate.R \
			--input {input} \
			{params.wildcard_flags} \
			{params.other_flags} \
			--output {output} \
			--threads {threads} \
			--log {log}
		"""

	
# Set rule targets
integrate = expand(os.path.join(config["output_dir"], "{prefix}", "batch:{batch}_normalisation:{normalisation_method}_integration:{integration_method}", f"integrated.{config.get('format', 'qs')}"), prefix=[os.path.join(os.path.basename(os.path.dirname(merged)), os.path.basename(os.path.dirname(normalised))) for merged, normalised in itertools.product(merge, normalise)], batch="+".join(config.get("batch", ["orig.ident"])), normalisation_method=config.get("normalisation-method", "LogNormalize"), integration_method=config.get("integration-method", None))