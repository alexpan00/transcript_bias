# This rule creates a kallisto index from the transcriptome FASTA file for short-read quantification.
rule kallisto_index_sr:
    input:
        fasta=rules.transcriptome_fasta.output.output,
    output:
        output_indx=os.path.join(config["output_dir"], "kallisto_sr", "01_index", "transcripts.idx"),
    conda:
        ENVS + "/kallisto.yaml"
    threads: 8
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/kallisto_sr/index.log"
    benchmark:
        BENCHMARKS + "/kallisto_sr/index.txt"
    params:
        tmp_dir=os.path.join(config["output_dir"], "kallisto_sr", "01_index", "tmp")
    shell:
        "kallisto index --tmp {params.tmp_dir} -i {output.output_indx} -t {threads} {input.fasta} > {log} 2>&1"

# This rule runs kallisto quant to quantify transcript abundances from short-read data.
rule kallisto_quant_sr:
    input:
        index=rules.kallisto_index_sr.output.output_indx,
        fastq=lambda wildcards: list(SR_FASTQ.loc[wildcards.sample])
    output:
        abundance=os.path.join(config["output_dir"], "kallisto_sr", "02_quant", "{sample}","abundance.tsv"),
        output_dir=directory(os.path.join(config["output_dir"], "kallisto_sr", "02_quant", "{sample}"))
    conda:
        ENVS + "/kallisto.yaml"
    threads: 8
    resources:
        mem_mb=7000,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/kallisto_sr/{sample}/quant.log"
    benchmark:
        BENCHMARKS + "/kallisto_sr/{sample}/quant.txt"
    params:
        seed=config.get("seed", 42)
    shell:
        '''
        kallisto quant -t {threads} --seed {params.seed} -i {input.index} -o {output.output_dir} \
        --rf-stranded {input.fastq} > {log} 2>&1
        '''

# This rule merges the kallisto quantification results into a single count matrix.
rule kallisto_counts_sr:
    wildcard_constraints:
        experiment=config["experiment"]
    input:
        abundance=expand(os.path.join(config["output_dir"], "kallisto_sr", "02_quant", "{sample}","abundance.tsv"), sample=grouped_by_sample.keys()),
    output:
        counts=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "{experiment}_counts.tsv")
    conda:
        ENVS + "/kallisto.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{experiment}/kallisto_sr/counts.log"
    benchmark:
        BENCHMARKS + "/{experiment}/kallisto_sr/counts.txt"
    run:
        import pandas as pd
        import os
        try:
            counts = {}
            for abundance in input.abundance:
                sample = os.path.basename(os.path.dirname(abundance))
                df = pd.read_csv(abundance, sep="\t")
                counts[sample] = df.set_index("target_id")["est_counts"]
            count_matrix = pd.DataFrame(counts)
            count_matrix = count_matrix.fillna(0)
            count_matrix.to_csv(output.counts, sep="\t", index=True, index_label="transcript_id")
        except Exception:
            log_and_raise(log[0])

# This rule prepares the transcriptome stats for kallisto short-read analysis.
rule kallisto_sr_transcriptome_stats:
    input:
        transcriptome_stats=rules.transcriptome_stats.output.output,
        counts=rules.kallisto_counts_sr.output.counts,
        script=SCRIPTS + "/analysis/expressed_transcripts.R"
    output:
        transcriptome_stats=os.path.join(config["output_dir"], "kallisto_sr","NOIseq",
                                          "{experiment}_transcript_models.tsv"),
        structural_category=os.path.join(config["output_dir"], "kallisto_sr","NOIseq",
                                          "{experiment}_structural_category.tsv")
    log:
        LOGS + "/{experiment}/kallisto_sr/transcriptome_stats.log"
    benchmark:
        BENCHMARKS + "/{experiment}/kallisto_sr/transcriptome_stats.txt"
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    shell:
        '''
        ln -s -r {input.transcriptome_stats} {output.transcriptome_stats}
        Rscript {input.script} {input.counts} {output.structural_category} > {log} 2>&1
        '''

# This rule prepares the transcript to gene mapping file for kallisto short-read analysis.
# It creates a relative symbolic link to the mapping file generated in the preprocessing step,
# which is built from the reference annotation.
rule kallisto_transcript_to_gene_sr:
    input:
        t2g=rules.prepare_transcript_to_gene_map.output.transcript_to_gene
    output:
        t2g=os.path.join(config["output_dir"], "kallisto_sr", "NOIseq", "{experiment}_transcript_to_gene.tsv")
    conda:
        ENVS + "/gffread.yaml"
    threads: 1
    log:
        LOGS + "/{experiment}/kallisto_sr/t2g.log"
    benchmark:
        BENCHMARKS + "/{experiment}/kallisto_sr/t2g.txt"
    shell:
        "ln -s -r {input.t2g} {output.t2g} > {log} 2>&1"