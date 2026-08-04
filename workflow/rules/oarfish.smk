# This rule creates a transcriptome FASTA file from the reference annotation and genome.
rule transcriptome_fasta:
    input:
        annotation=config["reference_annotation"],
        genome=config["reference_genome"]
    output:
        output=os.path.join(config["output_dir"], "oarfish", "transcripts.fasta")
    conda:
        ENVS + "/gffread.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        f"{LOGS}/oarfish/{config["experiment"]}_transcriptome_fasta.log"
    benchmark:
        f"{BENCHMARKS}/oarfish/{config["experiment"]}_transcriptome_fasta.txt"
    shell:
        "gffread -w {output} -g {input.genome} {input.annotation} > {log} 2>&1"

# This rule calculates statistics for the transcriptome.
rule transcriptome_stats:
    input:
        fasta=rules.transcriptome_fasta.output.output
    output:
        output=os.path.join(config["output_dir"], "transcripts.tsv")
    conda:
        ENVS + "/seqkit.yaml"
    threads: 8
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        f"{LOGS}/oarfish/{config["experiment"]}_transcriptome_stats.log"
    benchmark:
        f"{BENCHMARKS}/oarfish/{config["experiment"]}_transcriptome_stats.txt"
    shell:
        "seqkit fx2tab -i -n -l -g -H {input.fasta} > {output} 2> {log}"

# This rule creates a minimap2 index for the transcriptome.
rule index_trasncriptome:
    input:
        fasta=rules.transcriptome_fasta.output.output
    output:
        output=os.path.join(config["output_dir"], "oarfish", "transcripts.mmi")
    conda:
        ENVS + "/align.yaml"
    threads: 1
    resources:
        mem_mb=index_memory,
        slurm_extra="'--qos=short'"
    log:
        f"{LOGS}/oarfish/{config["experiment"]}_index_transcriptome.log"
    benchmark:
        f"{BENCHMARKS}/oarfish/{config["experiment"]}_index_transcriptome.txt"
    params:
        preset=f'{"map-hifi" if config["data_type"] == "pacbio" else "map-ont"}'
    shell:
        "minimap2 -x {params} -d {output} {input.fasta} > {log} 2>&1"

# This rule aligns the long reads to the transcriptome using minimap2.
rule align_transcriptome:
    input:
        fastq=lambda wildcards: grouped_by_sample[wildcards.sample]["fastq"],
        transcriptome=rules.index_trasncriptome.output.output
    output:
        bam=os.path.join(config["output_dir"], "oarfish", "01_aligments","{sample}.bam")
    conda:
        ENVS + "/align.yaml"
    threads: 8
    resources:
        mem_mb=align_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/oarfish/" + config["experiment"] + "/{sample}_align_transcriptome.log"
    benchmark:
        BENCHMARKS + "/oarfish/" + config["experiment"] + "/{sample}_align_transcriptome.txt"
    params:
        preset=f'{"map-hifi" if config["data_type"] == "pacbio" else "map-ont"}'
    shell:
        "(minimap2 -t {threads} -ax {params} {input.transcriptome} {input.fastq} | samtools view -b -@4 -o {output}) 2> {log}"

# This rule quantifies transcript abundances using oarfish.
rule oarfish_quantify:
    input:
        bam=rules.align_transcriptome.output.bam
    output:
        output=os.path.join(config["output_dir"], "oarfish", "02_quantify", "{sample}.quant")
    conda:
        ENVS + "/oarfish.yaml"
    threads: 8
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/oarfish/" + config["experiment"] + "/{sample}_oarfish_quantify.log"
    benchmark:
        BENCHMARKS + "/oarfish/" + config["experiment"] + "/{sample}_oarfish_quantify.txt"
    params:
        prefix=os.path.join(config["output_dir"], "oarfish", "02_quantify", "{sample}")
    shell:
        "oarfish -j {threads} -a {input.bam} -o {params.prefix} --filter-group no-filters --model-coverage > {log} 2>&1"

# This rule merges the oarfish quantification results into a single count matrix.
rule oarfish_merge:
    input:
        quant=expand(os.path.join(config["output_dir"], "oarfish", "02_quantify", "{sample}.quant"), sample=metadata["sample"]),
        script=SCRIPTS + "/quantification/merge_oarfish.py"
    output:
        merged_quant=os.path.join(config["output_dir"], "oarfish", "NOIseq", "{experiment}_counts.tsv")
    conda:
        ENVS + "/kallisto_counts.yaml"
    threads: 8
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/oarfish/{experiment}_merge.log"
    benchmark:
        BENCHMARKS + "/oarfish/{experiment}_merge.txt"
    shell:
        "python3 {input.script} {output} {input.quant} > {log} 2>&1"

# This rule prepares the transcriptome stats for oarfish.
rule oarfish_transcriptome_stats:
    input:
        transcriptome_stats=rules.transcriptome_stats.output.output,
        counts=rules.oarfish_merge.output.merged_quant,
        script=SCRIPTS + "/analysis/expressed_transcripts.R"
    output:
        transcriptome_stats=os.path.join(config["output_dir"], "oarfish","NOIseq",
                                          "{experiment}_transcript_models.tsv"),
        structural_category=os.path.join(config["output_dir"], "oarfish","NOIseq",
                                          "{experiment}_structural_category.tsv")
    conda:
        ENVS + "/NOIseq.yaml"
    log:
        LOGS + "/oarfish/{experiment}_transcriptome_stats.log"
    benchmark:
        BENCHMARKS + "/oarfish/{experiment}_transcriptome_stats.txt"
    threads: 1
    shell:
        '''
        ln -s -r {input.transcriptome_stats} {output.transcriptome_stats}
        Rscript {input.script} {input.counts} {output.structural_category} > {log} 2>&1
        '''

# This rule prepares the transcript to gene mapping file for oarfish.
# It creates a relative symbolic link to the mapping file generated in the preprocessing step,
# which is built from the reference annotation.
rule oarfish_transcript_to_gene:
    input:
        t2g=rules.prepare_transcript_to_gene_map.output.transcript_to_gene
    output:
        t2g=os.path.join(config["output_dir"], "oarfish", "NOIseq", "{experiment}_transcript_to_gene.tsv")
    conda:
        ENVS + "/gffread.yaml"
    threads: 1
    log:
        LOGS + "/oarfish/{experiment}_t2g.log"
    benchmark:
        BENCHMARKS + "/oarfish/{experiment}_t2g.txt"
    shell:
        "ln -s -r {input.t2g} {output.t2g} > {log} 2>&1"