# Helper function to generate dynamic report outputs for NOIseq plots
def get_noiseq_report_outputs(subdir_template, category, factors, norm_method, plot_set="all"):
    outputs = {}
    
    all_plot_types = {
        # name: (filename_part, caption_file, subcategory, extra_labels)
        "lenplot": ("length", "length.rst", "length", {"bins": "yes", "hex": "no", "normalization": norm_method}),
        "gcplot": ("GC", "GC.rst", "GC", {"bins": "yes", "hex": "no", "normalization": norm_method}),
        #"lenplot_raw_point": ("length_point", "length.rst", "length", {"bins": "no", "hex": "no", "normalization": norm_method}),
        #"lenplot_raw_hex": ("length_hex", "length.rst", "length", {"bins": "no", "hex": "yes", "normalization": norm_method}),
        #"gcplot_raw_point": ("GC_point", "GC.rst", "GC", {"bins": "no", "hex": "no", "normalization": norm_method}),
        #"gcplot_raw_hex": ("GC_hex", "GC.rst", "GC", {"bins": "no", "hex": "yes", "normalization": norm_method}),
        #"lenplot_bins_boxplot": ("Length_boxplot", "length.rst", "length", {"bins": "yes", "hex": "boxplot", "normalization": norm_method}),
        #"lenplot_ridge": ("Length_ridge", "length.rst", "length", {"bins": "yes", "hex": "ridge", "normalization": norm_method}),
        #"gcplot_boxplot": ("GC_boxplot", "GC.rst", "GC", {"bins": "yes", "hex": "boxplot", "normalization": norm_method}),
        #"gcplot_ridge": ("GC_ridge", "GC.rst", "GC", {"bins": "yes", "hex": "ridge", "normalization": norm_method}),
    }

    if plot_set == "sr":
        sr_plots = ["lenplot", "gcplot"]
        plot_types = {k: all_plot_types[k] for k in sr_plots}
    else:
        plot_types = all_plot_types

    for name, (pattern, caption_file, subcat, extra_labels) in plot_types.items():
        for factor in factors:
            labels = {"factor": factor}
            labels.update(extra_labels)
            
            path_pattern = "{experiment}_" + pattern + "_" + factor + ".png"
            path = os.path.join(config["output_dir"], subdir_template, "NOIseq", norm_method, path_pattern)

            outputs[f"{name}_{factor}"] = report(
                path,
                caption=os.path.join(REPORT, caption_file),
                category=category,
                subcategory=subcat,
                labels=labels
            )
            
    return outputs


def get_sirv_contrast_reports():
    outputs = {}
    # load the factors file to get the contrasts
    factor_df = pd.read_csv(config["factors"])

    # get the second column as a set
    factor_column = sorted(list(set(factor_df[factor_df.columns[1]])))
    contrasts = []

    for i in range(len(factor_column)):
        for j in range(i+1, len(factor_column)):
            contrasts.append(f"{factor_column[i]}_{factor_column[j]}")
    for contrast in contrasts:
        path = os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_sirvs_" + f"{contrast}.png")
        outputs[f"{contrast}"] = report(
            path,
            caption=REPORT + "/sirv_contrast.rst",
            category="{tool}",
            subcategory="validation",
            labels={
                "type": "SIRVs",
                "all": contrast,
                "plot": "contrast",
                "normalization": "{normalization_method}"
            }
        )
    return outputs
# This rule performs differential expression analysis using NOISeq.
rule NOIseq:
    input:
        NOISeq_object=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}","{experiment}_NOIseq.rds"),
        script=SCRIPTS + "/NOIseq_analysis.R",
    output:
        check=os.path.join(config["output_dir"], "{tool}", "NOIseq", ".{experiment}_{normalization_method}"),
        corplot=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}","{experiment}_heatmap.png"),
             caption=REPORT + "/correlation.rst",
             category="{tool}",
             subcategory="correlation",
             labels={
                "experiment": "{experiment}",
                "figure": "{experiment}_heatmap.png",
                "normalization": "{normalization_method}"
            }
        ),
        pcaplot=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}","{experiment}_pca.png"),
             caption=REPORT + "/pca.rst",
             category="{tool}",
             subcategory="pca",
             labels={
                "experiment": "{experiment}",
                "figure": "{experiment} pca",
                "normalization": "{normalization_method}"
            }
        ),
        sqanti_relative=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}","{experiment}_SQ_relative.png"),
             caption=REPORT + "/sqanti.rst",
             category="{tool}",
             subcategory="SQANTI3",
             labels={
                "absolute": "FALSE",
                "type": "Structural categories",
                "normalization": "{normalization_method}"
            }
        ),
        sqanti_absolute=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}","{experiment}_SQ_total.png"),
             caption=REPORT + "/sqanti.rst",
             category="{tool}",
             subcategory="SQANTI3",
             labels={
                "absolute": "TRUE",
                "type": "Structural categories",
                "normalization": "{normalization_method}"
            }
        ),
        **get_noiseq_report_outputs("{tool}", "{tool}", config.get("report_factors", []), "{normalization_method}")
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{tool}/{experiment}_{normalization_method}_NOIseq.log"
    benchmark:
        BENCHMARKS + "/{tool}/{experiment}_{normalization_method}_NOIseq.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}"),
        factors=",".join(config.get("report_factors", [])),
        type="lr"
    shell:
        '''
        Rscript {input.script} {input.NOISeq_object} \
            {params.output_prefix} --factors {params.factors} --type {params.type} > {log} 2>&1 && touch {output.check}
        '''

rule NOIseq_object:
    input:
        counts=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{experiment}_counts.tsv"),
        factors=config["factors"],
        transcripts_stats=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{experiment}_transcript_models.tsv"),
        script=SCRIPTS + "/NOIseq_object.R",
        structural_category=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{experiment}_structural_category.tsv"),
        metadata_extended=os.path.join(config["output_dir"], "metadata_extended.tsv")
    output:
        check=os.path.join(config["output_dir"], "{tool}", "NOIseq", ".{experiment}"),
        noiseq_obj=os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{experiment}_NOIseq.rds"),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{tool}/{experiment}_NOIseq_object.log"
    benchmark:
        BENCHMARKS + "/{tool}/{experiment}_NOIseq_object.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{experiment}"),
        factors=",".join(config.get("report_factors", [])),
        min_count_condition=config.get("min_count_condition", 0)
    shell:
        '''
        Rscript {input.script} {input.counts} {input.factors} \
            {input.transcripts_stats} {input.structural_category} \
            {params.output_prefix} --factors {params.factors} \
            {input.metadata_extended} {params.min_count_condition} > {log} 2>&1 && touch {output.check}
        '''

rule NOIseq_object_sr:
    input:
        counts=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "{experiment}_counts.tsv"),
        factors=config["factors"],
        transcripts_stats=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "{experiment}_transcript_models.tsv"),
        script=SCRIPTS + "/NOIseq_object.R",
        structural_category=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "{experiment}_structural_category.tsv"),
        metadata_extended=os.path.join(config["output_dir"], "metadata_extended.tsv")
    output:
        check=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", ".{experiment}"),
        noiseq_obj=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "{experiment}_NOIseq.rds"),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/kallisto_sr/{experiment}_NOIseq_object.log"
    benchmark:
        BENCHMARKS + "/kallisto_sr/{experiment}_NOIseq_object.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "{experiment}"),
        factors=",".join(config.get("report_factors", [])),
        min_count_condition=config.get("min_count_condition", 0)
    shell:
        '''
        Rscript {input.script} {input.counts} {input.factors} \
            {input.transcripts_stats} {input.structural_category} \
            {params.output_prefix} --factors {params.factors} \
            {input.metadata_extended} {params.min_count_condition} > {log} 2>&1 && touch {output.check}
        '''

# This rule performs differential expression analysis on short-read data using NOIseq.
rule NOIseq_sr:
    input:
        NOISeq_object=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "TPM","{experiment}_NOIseq.rds"),
        script=SCRIPTS + "/NOIseq_analysis.R",
    output:
        check=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", ".sr_{experiment}"),
        corplot=report(
            os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "TPM", "{experiment}_heatmap.png"),
             caption=REPORT + "/correlation.rst",
             category="kallisto_sr",
             subcategory="correlation",
             labels={
                "experiment": "{experiment}",
                "figure": "{experiment}_heatmap.png"
            }
        ),
        pcaplot=report(
            os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "TPM", "{experiment}_pca.png"),
             caption=REPORT + "/pca.rst",
             category="kallisto_sr",
             subcategory="pca",
             labels={
                "experiment": "{experiment}",
                "figure": "{experiment} pca"
            }
        ),
        **get_noiseq_report_outputs("kallisto_sr", "kallisto_sr", config.get("report_factors", []), "TPM",plot_set="sr")
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/kallisto_sr/{experiment}_NOIseq.log"
    benchmark:
        BENCHMARKS + "/kallisto_sr/{experiment}_NOIseq.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "TPM", "{experiment}"),
        factors=",".join(config.get("report_factors", [])),
        type="sr"
    shell:
        '''
        Rscript {input.script} {input.NOISeq_object} \
            {params.output_prefix} --factors {params.factors} --type {params.type} > {log} 2>&1 && touch {output.check}
        '''

# This rule compares the long-read and short-read quantification results.
def get_short_read_input(wildcards):
    if wildcards.normalization_method == "ratio_counts":
        folder = "ratio_counts"
    else:
        folder = "TPM"
        
    return os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", folder, f"{wildcards.experiment}_NOIseq.rds")

rule long_vs_short:
    input:
        long=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}","{experiment}_NOIseq.rds"),
        short=get_short_read_input,
        script=SCRIPTS + "/long_vs_short.R"
    output:
        check=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", ".{experiment}_long_vs_short"),
        summary_sr=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_sr_summary.rds"),
        cor_all=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_cor_all.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="validation",
             labels={
                "type": "short-reads",
                "all": "yes",
                "plot": "correlation",
                "level": "transcript",
                "normalization": "{normalization_method}",
            }
        ),
        cor_common=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_cor_common.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="validation",
             labels={
                "type": "short-reads",
                "all": "no",
                "plot": "correlation",
                "level": "transcript",
                "normalization": "{normalization_method}",
            }
        ),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/long_vs_short/{tool}/{normalization_method}/{experiment}_long_vs_short.log"
    benchmark:
        BENCHMARKS + "/long_vs_short/{tool}/{normalization_method}/{experiment}_long_vs_short.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}"),
    shell:
        '''
        Rscript {input.script} {input.long} {input.short} \
            {params.output_prefix} > {log} 2>&1 && touch {output.check}
        '''

# This rule compares the long-read and short-read quantification results for 
# a subset of transcripts.
rule long_vs_short_tusco:
    input:
        long=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}","{experiment}_NOIseq.rds"),
        short=get_short_read_input,
        script=SCRIPTS + "/compute_subset_correlation.R",
        tusco_list=config.get("tusco_list", ""),
    output:
        summary_tusco=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_tusco_summary.rds"),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/long_vs_short/{tool}/{normalization_method}/{experiment}_long_vs_short_tusco.log"
    benchmark:
        BENCHMARKS + "/long_vs_short/{tool}/{normalization_method}/{experiment}_long_vs_short_tusco.txt"
    shell:
        '''
        Rscript {input.script} {input.long} {input.short} \
            {input.tusco_list} {output.summary_tusco} > {log} 2>&1
        '''

# This rule compares the long-read and short-read quantification results.
def get_short_read_gene_level_input(wildcards):
    suffix = "gene_level_NOIseq.rds"
    if wildcards.normalization_method == "ratio_counts":
        folder = "ratio_counts"
    else:
        folder = "TPM"
        suffix = "gene_level_expression.rds"
        
    return os.path.join(config["output_dir"], 
            "kallisto_sr", 
            "NOIseq", 
            folder, 
            f"{wildcards.experiment}_{suffix}")

def get_long_read_gene_level_input(wildcards):
    suffix = "gene_level_NOIseq.rds"
    if wildcards.normalization_method == "ratio_counts":
        folder = "ratio_counts"
    else:
        folder = "TPM"
        suffix = "gene_level_expression.rds"
        
    return os.path.join(
        config["output_dir"], 
        wildcards.tool,
        "NOIseq", 
        folder, 
        f"{wildcards.experiment}_{suffix}"
    )


rule long_vs_short_gene_level:
    input:
        long=get_long_read_gene_level_input,
        short=get_short_read_gene_level_input,
        script=SCRIPTS + "/long_vs_short.R"
    output:
        check=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", ".{experiment}_gene_level_long_vs_short"),
        summary_sr=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_gene_level_sr_summary.rds"),
        cor_all=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_gene_level_cor_all.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="validation",
             labels={
                "type": "short-reads",
                "all": "yes",
                "plot": "correlation",
                "level": "gene",
                "normalization": "{normalization_method}",
            }
        ),
        cor_common=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_gene_level_cor_common.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="validation",
             labels={
                "type": "short-reads",
                "all": "no",
                "plot": "correlation",
                "level": "gene",
                "normalization": "{normalization_method}",
            }
        ),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/long_vs_short/{tool}/{normalization_method}/{experiment}_gene_level_long_vs_short.log"
    benchmark:
        BENCHMARKS + "/long_vs_short/{tool}/{normalization_method}/{experiment}_gene_level_long_vs_short.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_gene_level"),
    shell:
        '''
        Rscript {input.script} {input.long} {input.short} \
            {params.output_prefix} > {log} 2>&1 && touch {output.check}
        '''

# This rule validates the long-read quantification results using SIRV spike-ins.
def get_sirv_input(wildcards):
    suffix = "SIRV_NOIseq.rds"
    if wildcards.normalization_method == "ratio_correction":
        folder = "ratio_correction"
    elif wildcards.normalization_method == "ratio_counts":
        folder = "ratio_counts"
    else:
        suffix = "SIRV_counts.rds"
        folder = ""

    return os.path.join(config["output_dir"], "SIRVs", "NOIseq", folder, f"{wildcards.experiment}_{suffix}")

rule SIRV_validation:
    input:
        long=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}","{experiment}_NOIseq.rds"),
        sirv=lambda wildcards: get_sirv_input,
        script=SCRIPTS + "/SIRV_analysis.R",
        sirvs_info=config.get("sirvs_info", "sirvs_info.csv"),
    output:
        check=os.path.join(config["output_dir"], "{tool}", "NOIseq", ".{experiment}_{normalization_method}_SIRV"),
        summary_sirv=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_sirv_summary.rds"),
        cor_common=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_sirv_cor_common.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="validation",
             labels={
                "type": "SIRVs",
                "all": "no",
                "plot": "correlation",
                "normalization": "{normalization_method}"
            }
        ),
        sirv_box=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_sirv_boxplot.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="validation",
             labels={
                "type": "SIRVs",
                "all": "no",
                "plot": "boxplot",
                "normalization": "{normalization_method}"
            }
        ),
        sirv_len=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_sirv_vs_len.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="validation",
             labels={
                "type": "SIRVs",
                "all": "no",
                "plot": "lengthbias",
                "normalization": "{normalization_method}"
            }
        ),
        sirv_len_gt=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_sirv_vs_len_gt.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="validation",
             labels={
                "type": "SIRVs",
                "all": "no",
                "plot": "lengthbias_gt",
                "normalization": "{normalization_method}"
            }
        ),
        sirv_detection=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_sirv_detection.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="validation",
             labels={
                "type": "SIRVs",
                "all": "yes",
                "plot": "detection",
                "normalization": "{normalization_method}"
            }
        ),
        sirv_detection_10_counts=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_sirv_detection_10_reads.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="validation",
             labels={
                "type": "SIRVs",
                "all": "yes",
                "plot": "detection_10_counts",
                "normalization": "{normalization_method}"
            }
        ),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/long_vs_short/{tool}/{normalization_method}/{experiment}_long_vs_short_sirv.log"
    benchmark:
        BENCHMARKS + "/long_vs_short/{tool}/{normalization_method}/{experiment}_long_vs_short_sirv.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}"),
    shell:
        '''
        Rscript {input.script} {input.long} {input.sirv} {input.sirvs_info} \
            {params.output_prefix} > {log} 2>&1 && touch {output.check}
        '''


# This rule validates the long-read quantification results using ERCC spike-ins.
def get_ercc_input(wildcards):
    suffix = "ERCC_NOIseq.rds"
    if wildcards.normalization_method == "ratio_correction":
        folder = "ratio_correction"
    elif wildcards.normalization_method == "ratio_counts":
        folder = "ratio_counts"
    else:
        suffix = "ERCC_counts.rds"
        folder = ""

    return os.path.join(config["output_dir"], "SIRVs", "NOIseq", folder, f"{wildcards.experiment}_{suffix}")


rule ERCC_validation:
    input:
        long=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}","{experiment}_NOIseq.rds"),
        ercc=lambda wildcards: get_ercc_input,
        script=SCRIPTS + "/ERCC_analysis.R",
    output:
        summary_ercc=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_ercc_summary.rds"),
        cor_common=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_ercc_cor_common.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="validation",
             labels={
                "type": "ERCCs",
                "all": "no",
                "plot": "correlation",
                "normalization": "{normalization_method}"
            }
        ),
        ercc_detection=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_ercc_detection.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="validation",
             labels={
                "type": "ERCCs",
                "all": "yes",
                "plot": "detection",
                "normalization": "{normalization_method}"
            }
        ),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/long_vs_short/{tool}/{normalization_method}/{experiment}_long_vs_short_ercc.log"
    benchmark:
        BENCHMARKS + "/long_vs_short/{tool}/{normalization_method}/{experiment}_long_vs_short_ercc.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}"),
    shell:
        '''
        Rscript {input.script} {input.long} {input.ercc} \
            {params.output_prefix} > {log} 2>&1
        '''

# This rule creates a SQANTI catgories plot per sample. It also creates a csv
# file to later do the plot of the the different tools. It only makes sense to
# run this rule with the long-read data without normalization, since the SQANTI
# categories are not affected by normalization.
rule SQANTI_analysis:
    input:
        long=os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw","{experiment}_NOIseq.rds"),
        script=SCRIPTS + "/SQANTI_analysis.R",
    output:
        summary=os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{experiment}_SQ_counts_per_sample.csv"),
        summary_global=os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{experiment}_SQ_counts_global.csv"),
        sqanti_sample_pct=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw","{experiment}_SQ_stacked_pct.png"),
             caption=REPORT + "/sqanti.rst",
             category="{tool}",
             subcategory="SQANTI3",
             labels={
                "absolute": "False",
                "type": "Sample",
                "normalization": "raw"
            }
        ),
        sqanti_sample_cnt=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw","{experiment}_SQ_stacked_cnt.png"),
             caption=REPORT + "/sqanti.rst",
             category="{tool}",
             subcategory="SQANTI3",
             labels={
                "absolute": "True",
                "type": "Sample",
                "normalization": "raw"
            }
        ),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/sqanti/{tool}/{experiment}_sqanti_sample.log"
    benchmark:
        BENCHMARKS + "/sqanti/{tool}/{experiment}_sqanti_sample.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{experiment}"),
    shell:
        '''
        Rscript {input.script} {input.long} \
            {params.output_prefix} > {log} 2>&1
        '''

# This rule performs contrasts between the different submixes of SIRVs.
# It takes the SIRVs info and the long-reads quantification in the form of a NOIseq object
# and produces a plot for each contrast with the expected FC vs the observed FC.
rule SIRV_contrasts:
    input:
        long=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}","{experiment}_NOIseq.rds"),
        script=SCRIPTS + "/SIRV_contrasts.R",
        sirvs_info=config.get("sirvs_info", "sirvs_info.csv"),
    output:
        check=os.path.join(config["output_dir"], "{tool}", "NOIseq", ".{experiment}_{normalization_method}_SIRV_contrasts"),
        **get_sirv_contrast_reports(),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/sirvs/{tool}/{experiment}_{normalization_method}_SIRV_contrasts.log"
    benchmark:
        BENCHMARKS + "/sirvs/{tool}/{experiment}_{normalization_method}_SIRV_contrasts.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_sirvs"),
    shell:
        '''
        Rscript {input.script} {input.long} {input.sirvs_info} \
            {params.output_prefix} > {log} 2>&1 && touch {output.check}
       '''

# compute the the sensitivity in the detection of the SIRVs
rule sirv_sensitivity:
    input:
        long=os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{experiment}_NOIseq.rds"),
        sirv_object=os.path.join(config["output_dir"], "SIRVs", "NOIseq", "{experiment}_SIRV_counts.rds"),
        script=SCRIPTS + "/sirv_sensitivity.R",
    output:
        tool_sirv_sensitivity=os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{experiment}_sirv_sensitivity.csv"),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/sirvs/{tool}/{experiment}_SIRV_sensitivity.log"
    benchmark:
        BENCHMARKS + "/sirvs/{tool}/{experiment}_SIRV_sensitivity.txt"
    shell:
        '''
        Rscript {input.script} {input.long} {input.sirv_object} \
            {output.tool_sirv_sensitivity} > {log} 2>&1
        '''

# compute the the sensitivity in the detection of the ERCCs
rule ercc_sensitivity:
    input:
        long=os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{experiment}_NOIseq.rds"),
        ercc_object=os.path.join(config["output_dir"], "SIRVs", "NOIseq", "{experiment}_ERCC_counts.rds"),
        script=SCRIPTS + "/sirv_sensitivity.R",
    output:
        tool_ercc_sensitivity=os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{experiment}_ercc_sensitivity.csv"),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/sirvs/{tool}/{experiment}_ERCC_sensitivity.log"
    benchmark:
        BENCHMARKS + "/sirvs/{tool}/{experiment}_ERCC_sensitivity.txt"
    shell:
        '''
        Rscript {input.script} {input.long} {input.ercc_object} \
            {output.tool_ercc_sensitivity} > {log} 2>&1
       '''

# This rule explores the relationship between the number of isoforms per gene and gene expression.
rule isoforms_per_gene:
    input:
        long=rules.NOIseq_object.output.noiseq_obj,
        classification=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{experiment}_classification.txt"),
        script=SCRIPTS + "/expression_vs_n_isoforms.R",
    output:
        check=os.path.join(config["output_dir"], "{tool}", "NOIseq", ".{experiment}_iso_per_gene"),
        gene_expression_vs_isoforms=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{experiment}_gene_expression_vs_n_isoforms.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="random ideas",
             labels={
                "plot": "total_expression_vs_n_isoforms",
            }
        ),
        mean_expression_vs_isoforms=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{experiment}_mean_isoform_expression_vs_n_isoforms.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="random ideas",
             labels={
                "plot": "mean_expression_vs_n_isoforms",
            }
        ),
        max_expression_vs_isoforms=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{experiment}_max_isoform_expression_vs_n_isoforms.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="random ideas",
             labels={
                "plot": "max_expression_vs_n_isoforms",
            }
        ),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/ideas/{tool}/{experiment}_long_vs_short_sirv.log"
    benchmark:
        BENCHMARKS + "/ideas/{tool}/{experiment}_long_vs_short_sirv.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{experiment}"),
    shell:
        '''
        Rscript {input.script} {input.long} {input.classification} \
            {params.output_prefix} > {log} 2>&1 && touch {output.check}
        '''

# This rule assesses the replicability of the experiment.
rule replicability:
    input:
        long=rules.NOIseq_object.output.noiseq_obj,
        script=SCRIPTS + "/replicability.R",
    output:
        check=os.path.join(config["output_dir"], "{tool}", "NOIseq", ".{experiment}_replicability"),
        replicability_expression=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{experiment}_replicability_Expression.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="Replicability",
             labels={
                "plot": "Expression",
            }
        ),
        replicability_gc=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{experiment}_replicability_GC.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="Replicability",
             labels={
                "plot": "GC",
            }
        ),
        replicability_length=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{experiment}_replicability_Length.png"),
             caption=REPORT + "/long_vs_short.rst",
             category="{tool}",
             subcategory="Replicability",
             labels={
                "plot": "Length",
            }
        ),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/replicability/{tool}/{experiment}_long_vs_short_sirv.log"
    benchmark:
        BENCHMARKS + "/replicability/{tool}/{experiment}_long_vs_short_sirv.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{experiment}"),
    shell:
        '''
        Rscript {input.script} {input.long} \
            {params.output_prefix} > {log} 2>&1 && touch {output.check}
        '''

# This rule generates coverage plots.
rule coverage:
    input:
        bam=rules.align_transcriptome.output.bam,
        script=SCRIPTS + "/coverage.R",
        gclen=rules.transcriptome_stats.output.output,
    output:
        coverage_plot=report(
            os.path.join(config["output_dir"], "coverage", "{sample}_coverage.png"),
             caption=REPORT + "/coverage.rst",
             category="Reads",
             subcategory="Coverage",
             labels={
                "sample": "{sample}",
            }
        ),
    conda:
        ENVS + "/coverage.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/coverage/{sample}_coverage.log"
    benchmark:
        BENCHMARKS + "/coverage/{sample}_coverage.txt"
    shell:
        '''
        Rscript {input.script} {input.bam} {input.gclen} \
            {output.coverage_plot} > {log} 2>&1
        '''

rule normalization:
    input:
        noiseq_obj=rules.NOIseq_object.output.noiseq_obj,
        script=SCRIPTS + "/normalization.R",
        length_normalization_params=lambda wildcards: [config["length_normalization_params"]] if wildcards.norm_method == "read_density" else [],
    output:
        norm_counts = os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{experiment}_NOIseq.rds"),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{tool}/{experiment}_{norm_method}_normalization.log"
    benchmark:
        BENCHMARKS + "/{tool}/{experiment}_{norm_method}_normalization.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{experiment}"),
        post_filtering=config.get("post_filtering", 0),
    shell:
        '''
        Rscript {input.script} {input.noiseq_obj} \
            {params.output_prefix} {wildcards.norm_method} \
            {params.post_filtering} {input.length_normalization_params} > {log} 2>&1
        '''

# sr gene level expression normalization
rule short_read_gene_level_tpm:
    input:
        noiseq_obj=rules.NOIseq_object_sr.output.noiseq_obj,
        tx2gene_file=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "{experiment}_transcript_to_gene.tsv"),
        script_norm=SCRIPTS + "/normalization.R",
        script_gene_level_expr=SCRIPTS + "/get_gene_level_expression.R",
    output:
        norm_counts = os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "TPM", "{experiment}_NOIseq.rds"),
        gene_level_norm_counts = os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "TPM", "{experiment}_gene_level_expression.rds"),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/kallisto_sr/{experiment}_TPM_normalization.log"
    benchmark:
        BENCHMARKS + "/kallisto_sr/{experiment}_TPM_normalization.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "TPM", "{experiment}"),
        post_filtering=config.get("post_filtering", 0),
    shell:
        '''
        # First do TPM normalization at the transcript level
        Rscript {input.script_norm} {input.noiseq_obj} \
            {params.output_prefix} "TPM" {params.post_filtering} > {log} 2>&1
        # Then compute the gene level expression
        Rscript {input.script_gene_level_expr} {output.norm_counts} {input.tx2gene_file} \
            {output.gene_level_norm_counts} >> {log} 2>&1
        '''

# get long read gene level expression for all the normalization methods
rule long_read_gene_level_expression:
    input:
        noiseq_objs=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{experiment}_NOIseq.rds"),
        tx2gene_file=os.path.join(config["output_dir"],"{tool}", "NOIseq", "{experiment}_transcript_to_gene.tsv"),
        script=SCRIPTS + "/get_gene_level_expression.R",
    output:
        gene_level_norm_counts = os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{experiment}_gene_level_expression.rds"),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{tool}/{experiment}/{norm_method}_gene_level_expression.log"
    benchmark:
        BENCHMARKS + "/{tool}/{experiment}/{norm_method}_gene_level_expression.txt"
    shell:
        '''
        Rscript {input.script} {input.noiseq_objs} {input.tx2gene_file} \
            {output.gene_level_norm_counts} > {log} 2>&1
        '''
# create a summary plot for the dataset showing the effect of length normalization
# use the style of longbench 3c figure https://www.biorxiv.org/content/biorxiv/early/2025/09/12/2025.09.11.675724/F3.large.jpg?width=800&height=600&carousel=1
rule summary_normalization_length:
    input:
        noiseq_objs=expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{{experiment}}_NOIseq.rds"), norm_method = normalization_methods, tool=user_tools),
        script=SCRIPTS + "/summary_length.R",
    output:
        summary_df = os.path.join(config["output_dir"], "{experiment}_normalization_length_summary.csv"),
        summary_length_plot = report(
            os.path.join(config["output_dir"], "{experiment}_normalization_length_sum.png"),
            caption=REPORT + "/summary_length.rst",
            category="Summary",
            subcategory="Length normalization",
            labels={
                "plot": "normalization_length_sum",
            }
        ),
        mean_length_plot = report(
            os.path.join(config["output_dir"], "{experiment}_normalization_length_mean.png"),
            caption=REPORT + "/summary_length.rst",
            category="Summary",
            subcategory="Length normalization",
            labels={
                "plot": "normalization_length_mean",
            }
        ),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{experiment}/normalization_summary_length.log"
    benchmark:
        BENCHMARKS + "/{experiment}_normalization_summary_length.txt"
    params:
        noiseq_objs=",".join(expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{{experiment}}_NOIseq.rds"), norm_method = normalization_methods, tool=user_tools)),
    shell:
        '''
        Rscript {input.script} {params.noiseq_objs} {output.summary_length_plot} \
            {output.mean_length_plot} {output.summary_df} > {log} 2>&1
        '''

rule generate_synthetic_mixtures:
    input:
        raw=os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{experiment}_NOIseq.rds"),
        mix_def=config.get("mixture_definition", "mixture_definition.csv"),
        script=SCRIPTS + "/generate_synthetic_mixtures.R"
    output:
        syn_raw=os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{experiment}_synthetic_NOIseq.rds")
    conda:
        ENVS + "/NOIseq.yaml"
    log:
        LOGS + "/mixture_generation/{tool}/{experiment}_mixture_generation.log"
    benchmark:
        BENCHMARKS + "/mixture_generation/{tool}/{experiment}_mixture_generation.txt"
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    shell:
        '''
        Rscript {input.script} {input.raw} {input.mix_def} \
            {output.syn_raw} > {log} 2>&1
        '''

rule normalize_synthetic_mixtures:
    input:
        syn_raw=rules.generate_synthetic_mixtures.output.syn_raw,
        script=SCRIPTS + "/normalization.R",
        length_normalization_params=lambda wildcards: [config["length_normalization_params"]] if wildcards.normalization_method == "read_density" else [],
    output:
        syn_norm=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_synthetic_NOIseq.rds")
    conda:
        ENVS + "/NOIseq.yaml"
    params:
        output_prefix=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_synthetic"),
        post_filtering=config.get("post_filtering", 0),
    log:
        LOGS + "/mixture_normalization/{tool}/{normalization_method}/{experiment}_mixture_normalization.log"
    benchmark:
        BENCHMARKS + "/mixture_normalization/{tool}/{normalization_method}/{experiment}_mixture_normalization.txt"
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    shell:
        '''
        Rscript {input.script} {input.syn_raw} {params.output_prefix} \
            {wildcards.normalization_method} {params.post_filtering} \
            {input.length_normalization_params} > {log} 2>&1
        '''

rule mixture_correlation:
    input:
        obs_norm=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}","{experiment}_NOIseq.rds"),
        syn_norm=rules.normalize_synthetic_mixtures.output.syn_norm,
        mix_def=config.get("mixture_definition", "mixture_definition.csv"),
        script=SCRIPTS + "/mixture_correlation.R"
    output:
        check=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", ".{experiment}_mixture_correlation"),
        summary=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_mixture_summary.rds"),
        scatter=report(
            os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}_mixture_correlation.png"),
             caption=REPORT + "/mixture_correlation.rst",
             category="{tool}",
             subcategory="validation",
             labels={
                "type": "mixture",
                "plot": "correlation",
                "normalization": "{normalization_method}",
            }
        ),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/mixture_correlation/{tool}/{normalization_method}/{experiment}_mixture_correlation.log"
    benchmark:
        BENCHMARKS + "/mixture_correlation/{tool}/{normalization_method}/{experiment}_mixture_correlation.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{normalization_method}", "{experiment}"),
    shell:
        '''
        Rscript {input.script} {input.obs_norm} {input.syn_norm} {input.mix_def} \
            {params.output_prefix} > {log} 2>&1 && touch {output.check}
        '''

rule summary_normalization_mixtures:
    input:
        noiseq_objs=expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{{experiment}}_mixture_summary.rds"), norm_method = normalization_methods, tool=user_tools),
        script=SCRIPTS + "/summary_mixtures.R",
    output:
        summary_df = os.path.join(config["output_dir"], "{experiment}_normalization_mixture_summary.csv"),
        summary_correlation = report(
            os.path.join(config["output_dir"], "{experiment}_normalization_mixtures_correlation.png"),
            caption=REPORT + "/summary_mixture.rst",
            category="Summary",
            subcategory="Mixture normalization",
            labels={
                "plot": "normalization_correlation",
            }
        ),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{experiment}/normalization_summary_mixtures.log"
    benchmark:
        BENCHMARKS + "/{experiment}_normalization_summary_mixtures.txt"
    params:
        noiseq_objs=",".join(expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{{experiment}}_mixture_summary.rds"), norm_method = normalization_methods, tool=user_tools)),
    shell:
        '''
        Rscript {input.script} {params.noiseq_objs} {output.summary_correlation} \
            {output.summary_df} > {log} 2>&1
        '''

rule summary_normalization_sirv:
    input:
        noiseq_objs=expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{{experiment}}_sirv_summary.rds"), norm_method = normalization_methods, tool=user_tools),
        script=SCRIPTS + "/summary_sirv.R",
    output:
        summary_df = os.path.join(config["output_dir"], "{experiment}_normalization_sirv_summary.csv"),
        summary_correlation = report(
            os.path.join(config["output_dir"], "{experiment}_normalization_sirvs_correlation.png"),
            caption=REPORT + "/summary_sirv.rst",
            category="Summary",
            subcategory="SIRV normalization",
            labels={
                "plot": "normalization_correlation",
            }
        ),
        summary_error = report(
            os.path.join(config["output_dir"], "{experiment}_normalization_sirvs_error.png"),
            caption=REPORT + "/summary_sirv.rst",
            category="Summary",
            subcategory="SIRV normalization",
            labels={
                "plot": "normalization_error",
            }
        ),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{experiment}/normalization_summary_sirvs.log"
    benchmark:
        BENCHMARKS + "/{experiment}_normalization_summary_sirvs.txt"
    params:
        noiseq_objs=",".join(expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{{experiment}}_sirv_summary.rds"), norm_method = normalization_methods, tool=user_tools)),
    shell:
        '''
        Rscript {input.script} {params.noiseq_objs} {output.summary_correlation} \
            {output.summary_error} {output.summary_df} > {log} 2>&1
        '''

rule summary_sensitivity_sirv:
    input:
        sensitivity_objs=expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{{experiment}}_sirv_sensitivity.csv"), tool=user_tools),
        script=SCRIPTS + "/summary_sensitivity.R",
    output:
        summary_df = os.path.join(config["output_dir"], "{experiment}_sensitivity_sirv_summary.csv"),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{experiment}/sensitivity_summary_sirv.log"
    benchmark:
        BENCHMARKS + "/{experiment}/sensitivity_summary_sirv.txt"
    params:
        sensitivity_objs=",".join(expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{{experiment}}_sirv_sensitivity.csv"), tool=user_tools)),
    shell:
        '''
        Rscript {input.script} {params.sensitivity_objs} \
            {output.summary_df} > {log} 2>&1
        '''

rule summary_normalization_ercc:
    input:
        noiseq_objs=expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{{experiment}}_ercc_summary.rds"), norm_method = normalization_methods, tool=user_tools),
        script=SCRIPTS + "/summary_sirv.R",
    output:
        summary_df = os.path.join(config["output_dir"], "{experiment}_normalization_ercc_summary.csv"),
        summary_correlation = report(
            os.path.join(config["output_dir"], "{experiment}_normalization_ercc_correlation.png"),
            caption=REPORT + "/summary_sirv.rst",
            category="Summary",
            subcategory="ERCC normalization",
            labels={
                "plot": "normalization_correlation",
            }
        ),
        summary_error = report(
            os.path.join(config["output_dir"], "{experiment}_normalization_ercc_error.png"),
            caption=REPORT + "/summary_sirv.rst",
            category="Summary",
            subcategory="ERCC normalization",
            labels={
                "plot": "normalization_error",
            }
        ),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{experiment}/normalization_summary_ercc.log"
    benchmark:
        BENCHMARKS + "/{experiment}_normalization_summary_ercc.txt"
    params:
        noiseq_objs=",".join(expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{{experiment}}_ercc_summary.rds"), norm_method = normalization_methods, tool=user_tools)),
    shell:
        '''
        Rscript {input.script} {params.noiseq_objs} {output.summary_correlation} \
            {output.summary_error} {output.summary_df} > {log} 2>&1
        '''

rule summary_sensitivity_ercc:
    input:
        sensitivity_objs=expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{{experiment}}_ercc_sensitivity.csv"), tool=user_tools),
        script=SCRIPTS + "/summary_sensitivity.R",
    output:
        summary_df = os.path.join(config["output_dir"], "{experiment}_sensitivity_ercc_summary.csv"),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{experiment}/sensitivity_summary_ercc.log"
    benchmark:
        BENCHMARKS + "/{experiment}/sensitivity_summary_ercc.txt"
    params:
        sensitivity_objs=",".join(expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{{experiment}}_ercc_sensitivity.csv"), tool=user_tools)),
    shell:
        '''
        Rscript {input.script} {params.sensitivity_objs} \
            {output.summary_df} > {log} 2>&1
        '''

rule summary_normalization_sr:
    input:
        noiseq_objs=expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{{experiment}}_sr_summary.rds"), norm_method = normalization_methods, tool=user_tools),
        script=SCRIPTS + "/summary_sr.R",
    output:
        summary_df = os.path.join(config["output_dir"], "{experiment}_normalization_sr_summary.csv"),
        summary_correlation = report(
            os.path.join(config["output_dir"], "{experiment}_normalization_sr_correlation.png"),
            caption=REPORT + "/summary_sr.rst",
            category="Summary",
            subcategory="SR normalization",
            labels={
                "plot": "normalization_correlation",
                "level": "transcript",
            }
        ),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{experiment}/normalization_summary_sr.log"
    benchmark:
        BENCHMARKS + "/{experiment}_normalization_summary_sr.txt"
    params:
        noiseq_objs=",".join(expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{{experiment}}_sr_summary.rds"), norm_method = normalization_methods, tool=user_tools)),
        output_prefix=os.path.join(config["output_dir"], "{experiment}"),
    shell:
        '''
        Rscript {input.script} {params.noiseq_objs} {params.output_prefix} \
            {output.summary_df} > {log} 2>&1
        '''

rule summary_normalization_tusco:
    input:
        noiseq_objs=expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{{experiment}}_tusco_summary.rds"), norm_method = normalization_methods, tool=user_tools),
        script=SCRIPTS + "/summary_sr.R",
    output:
        summary_df = os.path.join(config["output_dir"], "{experiment}_normalization_tusco_summary.csv"),
        summary_correlation = report(
            os.path.join(config["output_dir"], "{experiment}_tusco_normalization_sr_correlation.png"),
            caption=REPORT + "/summary_sr.rst",
            category="Summary",
            subcategory="TUSCO normalization",
            labels={
                "plot": "normalization_correlation",
                "level": "transcript",
            }
        ),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{experiment}/normalization_summary_tusco.log"
    benchmark:
        BENCHMARKS + "/{experiment}_normalization_summary_tusco.txt"
    params:
        noiseq_objs=",".join(expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{{experiment}}_tusco_summary.rds"), norm_method = normalization_methods, tool=user_tools)),
        output_prefix=os.path.join(config["output_dir"], "{experiment}_tusco"),
    shell:
        '''
        Rscript {input.script} {params.noiseq_objs} {params.output_prefix} \
            {output.summary_df} > {log} 2>&1
        '''

rule summary_normalization_sr_gene_level:
    input:
        noiseq_objs=expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{{experiment}}_gene_level_sr_summary.rds"), norm_method = normalization_methods, tool=user_tools),
        script=SCRIPTS + "/summary_sr.R",
    output:
        summary_df = os.path.join(config["output_dir"], "{experiment}_gene_level_normalization_sr_summary.csv"),
        summary_correlation = report(
            os.path.join(config["output_dir"], "{experiment}_gene_level_normalization_sr_correlation.png"),
            caption=REPORT + "/summary_sr.rst",
            category="Summary",
            subcategory="SR normalization",
            labels={
                "plot": "normalization_correlation",
                "level": "gene",
            }
        ),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{experiment}/normalization_summary_sr_gene_level.log"
    benchmark:
        BENCHMARKS + "/{experiment}_normalization_summary_sr_gene_level.txt"
    params:
        noiseq_objs=",".join(expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "{norm_method}", "{{experiment}}_gene_level_sr_summary.rds"), norm_method = normalization_methods, tool=user_tools)),
        output_prefix=os.path.join(config["output_dir"], "{experiment}_gene_level"),
    shell:
        '''
        Rscript {input.script} {params.noiseq_objs} {params.output_prefix} \
            {output.summary_df} > {log} 2>&1
        '''
# Create a combined correlation overview plot (panel f) showing SR (transcript),
# SR (gene), mixture (optional) and SIRV (optional) correlations for every
# tool × normalization combination.  Mean correlations are computed per condition
# using Fisher z-transformation before averaging across samples.
rule summary_combined_correlation:
    input:
        sr_transcript_objs=rules.summary_normalization_sr.output.summary_df,
        sr_gene_objs=rules.summary_normalization_sr_gene_level.output.summary_df,
        mixture_objs=rules.summary_normalization_mixtures.output.summary_df if config.get("mixture_definition", False) else [],
        sirv_objs=rules.summary_normalization_sirv.output.summary_df if config.get("sirv_analysis", False) else [],
        ercc_objs=rules.summary_normalization_ercc.output.summary_df if config.get("ERCC_counts", False) else [],
        tusco_objs=rules.summary_normalization_tusco.output.summary_df if config.get("tusco_list", False) else [],
        script=SCRIPTS + "/summary_combined_correlation.R",
    output:
        summary_png=report(
            os.path.join(config["output_dir"], "{experiment}_combined_correlation.png"),
            caption=REPORT + "/summary_sr.rst",
            category="Summary",
            subcategory="Combined correlation",
            labels={
                "plot": "combined_correlation",
            }
        ),
        summary_csv=os.path.join(config["output_dir"], "{experiment}_combined_correlation_summary.csv"),
    conda:
        ENVS + "/combined_analysis.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{experiment}/summary_combined_correlation.log"
    benchmark:
        BENCHMARKS + "/{experiment}_summary_combined_correlation.txt"
    params:
        mixture_objs=lambda wildcards, input: input.mixture_objs if config.get("mixture_definition", False) else "none",
        sirv_objs=lambda wildcards, input: input.sirv_objs if config.get("sirv_analysis", False) else "none",
        ercc_objs=lambda wildcards, input: input.ercc_objs if config.get("ERCC_counts", False) else "none",
        tusco_objs=lambda wildcards, input: input.tusco_objs if config.get("tusco_list", False) else "none",
    shell:
        '''
        Rscript {input.script} {input.sr_transcript_objs} {input.sr_gene_objs} \
            {output.summary_png} {output.summary_csv} \
            {params.mixture_objs} {params.sirv_objs} {params.ercc_objs} {params.tusco_objs} > {log} 2>&1
        '''

# Create a summary plot for the SQANTI categories across tools
rule summary_sqanti:
    input:
        sqanti_csvs=expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{{experiment}}_SQ_counts_per_sample.csv"), tool=user_tools),
        script=SCRIPTS + "/summary_sqanti.R",
    output:
        summary_df = os.path.join(config["output_dir"], "{experiment}_sqanti_summary.csv"),
        summary_plot = report(
            os.path.join(config["output_dir"], "{experiment}_sqanti_summary.png"),
            caption=REPORT + "/sqanti.rst",
            category="Summary",
            subcategory="SQANTI3 summary",
            labels={
                "plot": "sqanti_summary_across_tools",
                "normalization": "raw"
            }
        ),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{experiment}/summary_sqanti_combined.log"
    benchmark:
        BENCHMARKS + "/{experiment}_summary_sqanti_combined.txt"
    params:
        sqanti_csvs=",".join(expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{{experiment}}_SQ_counts_per_sample.csv"), tool=user_tools)),
    shell:
        '''
        Rscript {input.script} {params.sqanti_csvs} {output.summary_plot} \
            {output.summary_df} > {log} 2>&1
        '''

# Create a summary plot for the SQANTI categories across tools but not by sample
# include the prefiltered and filtered trasncriptomes
rule summary_sqanti_global:
    input:
        sqanti_csvs=expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{{experiment}}_SQ_counts_global.csv"), tool=user_tools),
        categories_prefiltered=expand(os.path.join(config["output_dir"], "{tool}", "NOIseq","{{experiment}}_structural_category.tsv"), tool=user_tools),
        script=SCRIPTS + "/summary_sqanti_global.R",
    output:
        summary_df = os.path.join(config["output_dir"], "{experiment}_sqanti_summary_global.csv"),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{experiment}/summary_sqanti_global.log"
    benchmark:
        BENCHMARKS + "/{experiment}_summary_sqanti_global.txt"
    params:
        sqanti_csvs=",".join(expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{{experiment}}_SQ_counts_global.csv"), tool=user_tools)),
        categories_prefiltered=",".join(expand(os.path.join(config["output_dir"], "{tool}", "NOIseq","{{experiment}}_structural_category.tsv"), tool=user_tools)),
    shell:
        '''
        Rscript {input.script} {params.sqanti_csvs} {params.categories_prefiltered} \
            {output.summary_df} > {log} 2>&1
        '''
