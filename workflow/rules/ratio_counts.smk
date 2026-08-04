# get long read gene level expression for ratio_count normalization
rule long_read_gene_level_expression_ratio_counts:
    input:
        noiseq_objs=os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{experiment}_NOIseq.rds"),
        tx2gene_file=os.path.join(config["output_dir"],"{tool}", "NOIseq", "{experiment}_transcript_to_gene.tsv"),
        script=SCRIPTS + "/get_gene_level_expression.R",
        script_normalization=SCRIPTS + "/normalization.R",
    output:
        gene_level_norm_counts = os.path.join(config["output_dir"], "{tool}", "NOIseq", "raw", "{experiment}_gene_level_expression.rds"),
        gene_level_norm_counts_ratio = os.path.join(config["output_dir"], "{tool}", "NOIseq", "ratio_counts", "{experiment}_gene_level_NOIseq.rds"),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{tool}/{experiment}/ratio_counts_gene_level_expression.log"
    benchmark:
        BENCHMARKS + "/{tool}/{experiment}/ratio_counts_gene_level_expression.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "{tool}", "NOIseq", "ratio_counts", "{experiment}_gene_level"),
        post_filtering=config.get("post_filtering", 0),
    shell:
        '''
        Rscript {input.script} {input.noiseq_objs} {input.tx2gene_file} \
            {output.gene_level_norm_counts} > {log} 2>&1
        Rscript {input.script_normalization} {output.gene_level_norm_counts} \
            {params.output_prefix} ratio_counts \
            {params.post_filtering} >> {log} 2>&1
        '''
    
rule ratio_counts_normalization_sr:
    input:
        noiseq_obj=rules.NOIseq_object_sr.output.noiseq_obj,
        script=SCRIPTS + "/normalization.R",
    output:
        norm_counts = os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "ratio_counts", "{experiment}_NOIseq.rds"),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/kallisto_sr/{experiment}_ratio_counts_normalization.log"
    benchmark:
        BENCHMARKS + "/kallisto_sr/{experiment}_ratio_counts_normalization.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "ratio_counts", "{experiment}"),
        post_filtering=config.get("post_filtering", 0),
    shell:
        '''
        Rscript {input.script} {input.noiseq_obj} \
            {params.output_prefix} ratio_counts \
            {params.post_filtering} > {log} 2>&1
        '''



# get short read gene level expression for ratio_counts normalization
rule shortread_gene_level_expression_ratio_counts:
    input:
        noiseq_objs=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "{experiment}_NOIseq.rds"),
        tx2gene_file=os.path.join(config["output_dir"],"kallisto_sr", "NOIseq", "{experiment}_transcript_to_gene.tsv"),
        script=SCRIPTS + "/get_gene_level_expression.R",
        script_normalization=SCRIPTS + "/normalization.R",
    output:
        gene_level_norm_counts = os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "{experiment}_gene_level_expression.rds"),
        gene_level_norm_counts_ratio = os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "ratio_counts", "{experiment}_gene_level_NOIseq.rds"),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/kallisto_sr/{experiment}_ratio_counts_gene_level_expression.log"
    benchmark:
        BENCHMARKS + "/kallisto_sr/{experiment}_ratio_counts_gene_level_expression.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "ratio_counts", "{experiment}_gene_level"),
        post_filtering=config.get("post_filtering", 0),
    shell:
        '''
        Rscript {input.script} {input.noiseq_objs} {input.tx2gene_file} \
            {output.gene_level_norm_counts} > {log} 2>&1
        Rscript {input.script_normalization} {output.gene_level_norm_counts} \
            {params.output_prefix} ratio_counts \
            {params.post_filtering} >> {log} 2>&1
        '''

# get SIRV expression for ratio_counts normalization
rule SIRV_expression_ratio_counts:
    input:
        sirv_obj=os.path.join(config["output_dir"], "SIRVs", "NOIseq", "{experiment}_SIRV_counts.rds"),
        script_normalization=SCRIPTS + "/normalization.R",
    output:
        sirv_norm_counts_ratio = os.path.join(config["output_dir"], "SIRVs", "NOIseq", "ratio_counts", "{experiment}_SIRV_NOIseq.rds"),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/SIRVs/NOIseq/{experiment}_ratio_counts.log"
    benchmark:
        BENCHMARKS + "/SIRVs/NOIseq/{experiment}_ratio_counts.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "SIRVs", "NOIseq", "ratio_counts", "{experiment}_SIRV"),
        post_filtering=0,
    shell:
        '''
        Rscript {input.script_normalization} {input.sirv_obj} \
            {params.output_prefix} ratio_counts \
            {params.post_filtering} > {log} 2>&1
        '''

# get ERCC expression for ratio_counts normalization
rule ERCC_expression_ratio_counts:
    input:
        ercc_obj=os.path.join(config["output_dir"], "SIRVs", "NOIseq", "{experiment}_ERCC_counts.rds"),
        script_normalization=SCRIPTS + "/normalization.R",
    output:
        ercc_norm_counts_ratio = os.path.join(config["output_dir"], "SIRVs", "NOIseq", "ratio_counts", "{experiment}_ERCC_NOIseq.rds"),
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/ERCC/NOIseq/{experiment}_ratio_counts.log"
    benchmark:
        BENCHMARKS + "/ERCC/NOIseq/{experiment}_ratio_counts.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "SIRVs", "NOIseq", "ratio_counts", "{experiment}_ERCC"),
        post_filtering=0,
    shell:
        '''
        Rscript {input.script_normalization} {input.ercc_obj} \
            {params.output_prefix} ratio_counts \
            {params.post_filtering} > {log} 2>&1
        '''

# This rule performs differential expression analysis on short-read data using NOIseq.
rule NOIseq_sr_ratio_counts_analysis:
    input:
        NOISeq_object=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "ratio_counts","{experiment}_NOIseq.rds"),
        script=SCRIPTS + "/NOIseq_analysis.R",
    output:
        check=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", ".sr_{experiment}_ratio_counts"),
        corplot=report(
            os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "ratio_counts", "{experiment}_heatmap.png"),
             caption=REPORT + "/correlation.rst",
             category="kallisto_sr",
             subcategory="correlation",
             labels={
                "experiment": "{experiment}",
                "figure": "{experiment}_heatmap.png",
                "normalization": "ratio_counts"
            }
        ),
        pcaplot=report(
            os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "ratio_counts", "{experiment}_pca.png"),
             caption=REPORT + "/pca.rst",
             category="kallisto_sr",
             subcategory="pca",
             labels={
                "experiment": "{experiment}",
                "figure": "{experiment} pca",
                "normalization": "ratio_counts"
            }
        ),
        **get_noiseq_report_outputs("kallisto_sr", "kallisto_sr", config.get("report_factors", []), "ratio_counts",plot_set="sr")
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/kallisto_sr/{experiment}_ratio_counts_NOIseq.log"
    benchmark:
        BENCHMARKS + "/kallisto_sr/{experiment}_ratio_counts_NOIseq.txt"
    params:
        output_prefix=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "ratio_counts", "{experiment}"),
        factors=",".join(config.get("report_factors", [])),
        type="sr"
    shell:
        '''
        Rscript {input.script} {input.NOISeq_object} \
            {params.output_prefix} --factors {params.factors} --type {params.type} > {log} 2>&1 && touch {output.check}
        '''

