##########################################################################################
# Snakemake rule for merging multisample multimodal counts and metadata
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
##########################################################################################

scripts_dir = config.get("scripts_dir", "resources/scripts")

# Define rule
rule merge:
	output: os.path.join(config["output_dir"], "merge", f"{{merge}}.{config.get('format', 'qs')}"), 
	log: os.path.abspath("logs/merge/{merge}.log")
	threads: min(len(samples), 4)
	params:
		script_path = scripts_dir if os.path.isabs(scripts_dir) else os.path.join(workflow.basedir, scripts_dir),
		samples = ";".join(samples),
		hdf5 = ";".join([info[sample]["hdf5"] for sample in samples]),
		metadata = ";".join([info[sample]["metadata"] for sample in samples]),
		fragments_flag = f"--fragments '{';'.join([info[sample].get('fragments', '') for sample in samples])}'" if any(info[sample].get("fragments", None) for sample in samples) else "",
		optional_flags = get_optional_flags(gene_types=";".join(config.get("gene-types", list())), control=config.get("control", None), binarize=config.get("binarize", None), bpcells=config.get("bpcells", None), downsample=config.get("downsample", None), exclude_control=config.get("exclude-control", None)),
	conda: "single_cell_multi"
	# envmodules: "R-cbrg"
	message: "Merging counts and metadata"
	shell:
		"""
		cd {params.script_path} && \
		./merge.R \
			--samples '{params.samples}' \
			--hdf5 '{params.hdf5}' \
			--metadata '{params.metadata}' \
			{params.fragments_flag} \
			{params.optional_flags} \
			--output {output} \
			--threads {threads} \
			--log {log}
		"""

	
# Set rule targets
merge = expand(os.path.join(config["output_dir"], "merge", f"b:{{binarize}}_B:{{bpcells}}_d:{{downsample}}_e:{{exclude_control}}.{config.get('format', 'qs')}"), binarize=config.get("binarize", False), bpcells=config.get("bpcells", False), downsample=config.get("downsample", False), exclude_control=config.get("exclude-control", False))