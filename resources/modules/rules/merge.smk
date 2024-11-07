##########################################################################################
# Snakemake rule for merging multisample multimodal counts and metadata
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
##########################################################################################

scripts_dir = config.get("scripts_dir", "resources/scripts")

# Define rule
rule merge:
	output: os.path.join(config["output_dir"], "peak-method:{peak_method}_b:{binarize}_d:{downsample}", f"merged.{config.get('format', 'qs')}"), 
	log: os.path.abspath(os.path.join("logs", "peak-method:{peak_method}_b:{binarize}_d:{downsample}", "merge.log"))
	threads: min(floor(len(samples) / 5) + 1, 4)
	resources:
		mem = lambda wildcards, threads: f"{threads * 25}GiB"
	params:
		script_path = scripts_dir if os.path.isabs(scripts_dir) else os.path.join(workflow.basedir, scripts_dir),
		samples = ";".join(samples),
		hdf5 = ";".join([info[sample]["hdf5"] for sample in samples]),
		metadata = ";".join([info[sample]["metadata"] for sample in samples]),
		fragments_flag = f"--fragments '{';'.join([info[sample].get('fragments', '') for sample in samples])}'" if any(info[sample].get("fragments", None) for sample in samples) else "",
		summits_flag = f"--summits '{';'.join([info[sample].get('summits', '') for sample in samples])}'" if any(info[sample].get("summits", None) for sample in samples) else "",
		wildcard_flags = lambda wildcards: get_optional_flags(peak_method=wildcards.peak_method if wildcards.peak_method != "None" else None, binarize=str_to_bool(wildcards.binarize), downsample=str_to_bool(wildcards.downsample)),
		other_flags = get_optional_flags(rna_assay=config.get("rna-assay", None), atac_assay=config.get("atac-assay", None), adt_assay=config.get("adt-assay", None), gene_types=";".join(config.get("gene-types", list())), control=config.get("control", None), exclude_control=config.get("exclude-control", None)),
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
			{params.summits_flag} \
			{params.wildcard_flags} \
			{params.other_flags} \
			--output {output} \
			--threads {threads} \
			--log {log}
		"""

	
# Set rule targets
merge = expand(os.path.join(config["output_dir"], "peak-method:{peak_method}_b:{binarize}_d:{downsample}", f"merged.{config.get('format', 'qs')}"), peak_method = config.get("peak-method", None), binarize=config.get("binarize", False), downsample=config.get("downsample", False))