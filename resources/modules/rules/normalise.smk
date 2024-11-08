##########################################################################################
# Snakemake rule for normalising multisample multimodal counts and finding variable
# features
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
# Requires outputs from resources/rules/merge.smk
##########################################################################################

scripts_dir = config.get("scripts_dir", "resources/scripts")

# Define rule
rule normalise:
	input: os.path.join(config["output_dir"], "{prefix}", f"merged.{config.get('format', 'qs')}"),
	output: os.path.join(config["output_dir"], "{prefix}", "n-features:{n_features}_clr:{clr}_M:{regress_percent_mitochondrial}_C:{regress_cell_cycle}", f"normalised.{config.get('format', 'qs')}"), 
	log: os.path.abspath("logs/{prefix}/n-features:{n_features}_clr:{clr}_M:{regress_percent_mitochondrial}_C:{regress_cell_cycle}/normalise.log")
	threads: min(floor(len(samples) / 5) + 1, 4)
	resources:
		mem = lambda wildcards, threads: f"{threads * 100}GiB"
	params:
		script_path = scripts_dir if os.path.isabs(scripts_dir) else os.path.join(workflow.basedir, scripts_dir),
		wildcard_flags = lambda wildcards: get_optional_flags(n_features=wildcards.n_features, clr=wildcards.clr, regress_percent_mitochondrial=str_to_bool(wildcards.regress_percent_mitochondrial), regress_cell_cycle=str_to_bool(wildcards.regress_cell_cycle)),
		other_flags = get_optional_flags(rna_assay = config.get("rna-assay", None), atac_assay = config.get("atac-assay", None), adt_assay = config.get("adt-assay", None), exp=config.get("exp", None), rna_filter_features=config.get("rna-filter-features", None), atac_filter_features=config.get("atac-filter-features", None)),
	conda: "single_cell_multi"
	# envmodules: "R-cbrg"
	message: "Normalising counts and finding variable features"
	shell:
		"""
		cd {params.script_path} && \
		./normalise.R \
			--input {input} \
			{params.wildcard_flags} \
			{params.other_flags} \
			--output {output} \
			--threads {threads} \
			--log {log}
		"""

	
# Set rule targets
normalise = expand(os.path.join(config["output_dir"], "{prefix}", "n-features:{n_features}_clr:{clr}_M:{regress_percent_mitochondrial}_C:{regress_cell_cycle}", f"normalised.{config.get('format', 'qs')}"), prefix=[os.path.basename(os.path.dirname(merged)) for merged in merge], n_features=config.get("n-features", 3000), clr=config.get("clr", "seurat"), regress_percent_mitochondrial=config.get("regress-percent-mitochondrial", False), regress_cell_cycle=config.get("regress-cell-cycle", False))