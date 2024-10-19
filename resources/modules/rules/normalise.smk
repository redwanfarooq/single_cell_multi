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
	input: merge,
	output: os.path.join(config["output_dir"], "normalise", f"{{normalise}}.{config.get('format', 'qs')}"), 
	log: os.path.abspath("logs/normalise/{normalise}.log")
	threads: min(len(samples), 4)
	params:
		script_path = scripts_dir if os.path.isabs(scripts_dir) else os.path.join(workflow.basedir, scripts_dir),
		optional_flags = get_optional_flags(rna_assay = config.get("rna-assay", None), atac_assay = config.get("atac-assay", None), adt_assay = config.get("adt-assay", None), exp=config.get("exp", None), n_features=config.get("n-features", None), rna_filter_features=config.get("rna-filter-features", None), atac_filter_features=config.get("atac-filter-features", None), clr=config.get("clr", None), regress_percent_mitochondrial=config.get("regress-percent-mitochondrial", None), regress_cell_cycle=config.get("regress-cell-cycle", None)),
	conda: "single_cell_multi"
	# envmodules: "R-cbrg"
	message: "Normalising counts and finding variable features"
	shell:
		"""
		cd {params.script_path} && \
		./normalise.R \
			--input {input} \
			{params.optional_flags} \
			--output {output} \
			--threads {threads} \
			--log {log}
		"""

	
# Set rule targets
normalise = expand(os.path.join(config["output_dir"], "normalise", f"n-features:{{n_features}}_clr:{{clr}}_M:{{regress_percent_mitochondrial}}_C:{{regress_cell_cycle}}.{config.get('format', 'qs')}"), n_features=config.get("n-features", 3000), clr=config.get("clr", "seurat"), regress_percent_mitochondrial=config.get("regress-percent-mitochondrial", False), regress_cell_cycle=config.get("regress-cell-cycle", False))