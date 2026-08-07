# This rule creates a kallisto index from the reference genome and annotation.
rule kallisto_index:
    input:
        fasta=config["reference_genome"],
        annotation=config["reference_annotation"]
    output:
        output_indx=os.path.join(config["output_dir"], "kallisto", "01_index", "transcripts.idx"),
        output_t2g=os.path.join(config["output_dir"], "kallisto", "01_index", "transcripts.t2g"),
        fasta=os.path.join(config["output_dir"], "kallisto", "01_index", "transcripts.fasta")
    conda:
        ENVS + "/kallisto.yaml"
    threads: 8
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/" + config["experiment"] + "/kallisto" + "/index.log"
    benchmark:
        BENCHMARKS + "/" + config["experiment"] + "/kallisto" + "/index.txt"
    params:
        tmp_dir=os.path.join(config["output_dir"], "kallisto", "01_index", "tmp") 
    shell:
        "kb ref -k 63 --tmp {params.tmp_dir} -i {output.output_indx} -g {output.output_t2g} -f1 {output.fasta} {input.fasta} {input.annotation} > {log} 2>&1"

# This rule runs kallisto bus to generate BUS files from long-read data.
rule kallisto_bus:
    input:
        index=rules.kallisto_index.output.output_indx,
        fastq=list(metadata.fastq)
    output:
        output_bus=os.path.join(config["output_dir"], "kallisto", "02_bus", "output.bus"),
        output_trans=os.path.join(config["output_dir"], "kallisto", "02_bus", "transcripts.txt"),
        output_mat=os.path.join(config["output_dir"], "kallisto", "02_bus", "matrix.ec"),
        output_flens=os.path.join(config["output_dir"], "kallisto", "02_bus", "flens.txt"),
        output_dir=directory(os.path.join(config["output_dir"], "kallisto", "02_bus"))
    conda:
        ENVS + "/kallisto.yaml"
    threads: 8
    resources:
        mem_mb=7000,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/" + config["experiment"] + "/kallisto" + "/bus.log"
    benchmark:
        BENCHMARKS + "/" + config["experiment"] + "/kallisto" + "/bus.txt"
    params:
        "--long --threshold=0.8 -x bulk"
    shell:
        "kallisto bus {params} -t {threads} -i {input.index} -o {output.output_dir} {input.fastq} > {log} 2>&1"

# This rule sorts the BUS file.
rule bus_sort:
    input:
        bus=rules.kallisto_bus.output.output_bus
    output:
        output_bus=os.path.join(config["output_dir"], "kallisto", "02_bus", "output.sorted.bus")
    conda:
        ENVS + "/kallisto.yaml"
    threads: 8
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/" + config["experiment"] + "/kallisto" + "/bus_sort.log"
    benchmark:
        BENCHMARKS + "/" + config["experiment"] + "/kallisto" + "/bus_sort.txt"
    shell:
        "bustools sort -t {threads} -o {output.output_bus} {input.bus} > {log} 2>&1"

    
# This rule counts the reads in the BUS file.
rule bus_count:
    input:
        bus=rules.bus_sort.output.output_bus,
        t2g=rules.kallisto_index.output.output_t2g,
        ec=rules.kallisto_bus.output.output_mat,
        trans=rules.kallisto_bus.output.output_trans
    output:
        output_mat=os.path.join(config["output_dir"], "kallisto", "03_bus_count", "output.mtx"),
        output_mat_ec=os.path.join(config["output_dir"], "kallisto", "03_bus_count", "output.ec.txt"),
        output_dir=directory(os.path.join(config["output_dir"], "kallisto", "03_bus_count"))
    conda:
        ENVS + "/kallisto.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/" + config["experiment"] + "/kallisto" + "/bus_count.log"
    benchmark:
        BENCHMARKS + "/" + config["experiment"] + "/kallisto" + "/bus_count.txt"
    params:
        "--cm -m"
    shell:
        """
        bustools count {input.bus} -t {input.trans} -e {input.ec} -g {input.t2g} {params} -o {output.output_dir}/ > {log} 2>&1
        """

# This rule performs transcript-level quantification using kallisto quant-tcc.
rule kallisto_tcc:
    input:
        index=rules.kallisto_index.output.output_indx,
        flens=rules.kallisto_bus.output.output_flens,
        ec=rules.bus_count.output.output_mat_ec,
        mtx=rules.bus_count.output.output_mat
    output:
        output_mtx=os.path.join(config["output_dir"], "kallisto", "04_tcc", "matrix.abundance.mtx"),
        output_dir=directory(os.path.join(config["output_dir"], "kallisto", "04_tcc")),
        trans_ids=os.path.join(config["output_dir"], "kallisto", "04_tcc", "transcripts.txt")
    conda:
        ENVS + "/kallisto.yaml"
    threads: 8
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/" + config["experiment"] + "/kallisto" + "/tcc.log"
    benchmark:
        BENCHMARKS + "/" + config["experiment"] + "/kallisto" + "/tcc.txt"
    params:
        lambda wildcards, threads: f'--long -P {"PacBio" if config["data_type"] == "pacbio" else "ONT"}' # TODO should be an input of th pipeline
    shell:       
        "kallisto quant-tcc -t {threads}  {params} -f {input.flens} -i {input.index} "
        "-e {input.ec} -o {output.output_dir} {input.mtx} > {log} 2>&1"

# This rule converts the kallisto output to a count matrix.
rule kallisto_counts:
    input:
        counts=rules.kallisto_tcc.output.output_mtx,
        transcripts_ids=rules.kallisto_tcc.output.trans_ids,
        metadata_csv=config["metadata"],
        script=SCRIPTS + "/quantification/counts_kallisto.py"
    output:
        counts=os.path.join(config["output_dir"], "kallisto","NOIseq", "{experiment}_counts.tsv")
    conda:
        ENVS + "/kallisto_counts.yaml"
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{experiment}/kallisto/counts.log"
    benchmark:
        BENCHMARKS + "/{experiment}/kallisto/counts.txt"
    shell:
        "python3 {input.script} -c {input.counts} -t {input.transcripts_ids} -m {input.metadata_csv} --output {output} > {log} 2>&1"

# This rule prepares the transcriptome stats for kallisto.
rule kallisto_transcriptome_stats:
    input:
        transcriptome_stats=rules.transcriptome_stats.output.output,
        counts=rules.kallisto_counts.output.counts,
        script=SCRIPTS + "/analysis/expressed_transcripts.R"
    output:
        transcriptome_stats=os.path.join(config["output_dir"], "kallisto","NOIseq",
                                          "{experiment}_transcript_models.tsv"),
        structural_category=os.path.join(config["output_dir"], "kallisto","NOIseq",
                                          "{experiment}_structural_category.tsv")
    conda:
        ENVS + "/NOIseq.yaml"
    log:
        LOGS + "/{experiment}/kallisto/transcriptome_stats.log"
    benchmark:
        BENCHMARKS + "/{experiment}/kallisto/transcriptome_stats.txt"
    threads: 1
    shell:
        '''
        ln -s -r {input.transcriptome_stats} {output.transcriptome_stats}
        Rscript {input.script} {input.counts} {output.structural_category} > {log} 2>&1
        '''

# This rule prepares the transcript to gene mapping file for kallisto.
# It creates a relative symbolic link to the mapping file generated in the preprocessing step,
# which is built from the reference annotation.
rule kallisto_transcript_to_gene:
    input:
        t2g=rules.prepare_transcript_to_gene_map.output.transcript_to_gene
    output:
        t2g=os.path.join(config["output_dir"], "kallisto", "NOIseq", "{experiment}_transcript_to_gene.tsv")
    conda:
        ENVS + "/gffread.yaml"
    threads: 1
    log:
        LOGS + "/{experiment}/kallisto/t2g.log"
    benchmark:
        BENCHMARKS + "/{experiment}/kallisto/t2g.txt"
    shell:
        "ln -s -r {input.t2g} {output.t2g} > {log} 2>&1"