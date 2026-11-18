# This rule creates a gene database for IsoQuant.
rule isoquant_db:
    input:
        annotation=config["reference_annotation"],
    output:
        db=os.path.join(config["output_dir"], "isoquant", "{experiment}","isoquant.db")
    conda:
        ENVS + "/isoquant.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{experiment}/isoquant/isoquant_db.log"
    benchmark:
        BENCHMARKS + "/{experiment}/isoquant/isoquant_db.txt"
    script:
        SCRIPTS + "/preprocessing/isoquant_db.py"


# This rule runs IsoQuant to perform transcript reconstruction and quantification.
rule isoquant:
    input:
        bam=lambda wildcards: grouped[wildcards.condition]["aligned"],
        index=lambda wildcards: grouped[wildcards.condition]["bai"],
        db=rules.isoquant_db.output.db
    output:
        isoquant_gtf=os.path.join(config["output_dir"], "isoquant","{experiment}","{condition}","{condition}.transcript_models.gtf"),
        isoquant_counts=os.path.join(config["output_dir"], "isoquant","{experiment}","{condition}","{condition}.transcript_model_grouped_counts.tsv")
    conda:
        ENVS + "/isoquant.yaml"
    threads: get_rule_resource(config, "isoquant", "threads", 8)
    resources:
        mem_mb=get_rule_resource(config, "isoquant", "mem_mb", extra_memory),
        slurm_extra=get_rule_resource(config, "isoquant", "slurm_extra", "'--qos=short'")
    log:
        LOGS + "/isoquant/{experiment}/{condition}_isoquant.log"
    benchmark:
        BENCHMARKS + "/isoquant/{experiment}/{condition}_isoquant.txt"
    params:
        genome=f"--reference {config['reference_genome']}",
        annotation= lambda w,input: f"--genedb {input.db}",
        data_type=f"--data_type {config['data_type']}",
        output_dir="--output " + os.path.join(config["output_dir"], "isoquant", "{experiment}"),
        prefix="--prefix {condition}",
        # Name the count columns with the metadata sample IDs rather than the BAM
        # basenames, so they match the factors file however the BAMs are named.
        labels=lambda w: "--labels " + " ".join(grouped[w.condition]["sample"]),
        # Set only when the condition has a single BAM; see the rename below.
        single_sample=lambda w: grouped[w.condition]["sample"][0] if len(grouped[w.condition]["sample"]) == 1 else ""
    shell:
        """
        isoquant.py -t {threads} --bam {input.bam} {params.genome} {params.annotation} {params.data_type} {params.output_dir} {params.prefix} {params.labels} > {log} 2>&1

        # With a single BAM, IsoQuant writes <condition>.transcript_model_counts.tsv
        # with a "#feature_id<TAB>count" header instead of the grouped file. Rename it
        # and put the sample ID in the header, as in the multi-BAM case.
        if [ -n "{params.single_sample}" ]; then
            single_counts=$(dirname {output.isoquant_counts})/{wildcards.condition}.transcript_model_counts.tsv
            awk -v s="{params.single_sample}" 'BEGIN {{OFS="\\t"}} NR == 1 && /^#feature_id/ {{print "#feature_id", s; next}} {{print}}' "$single_counts" > {output.isoquant_counts}
            rm "$single_counts"
        fi
        """

# This rule converts the IsoQuant GTF output to BED format.
rule isoquant_bed:
    input:
        gtf=rules.isoquant.output.isoquant_gtf,
        gtf2bed=rules.prepare_tama.output.gtf2bed,
        script=SCRIPTS + "/preprocessing/gtf2bed.sh"
    output:
        bed=os.path.join(config["output_dir"], "isoquant","tama","{experiment}","{condition}.bed")
    conda:
        ENVS + "/tama.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/isoquant/{experiment}/{condition}_isoquant_bed.log"
    benchmark:
        BENCHMARKS + "/isoquant/{experiment}/{condition}_isoquant_bed.txt"
    shell:
        "bash {input.script} {input.gtf} {output} {input.gtf2bed} > {log} 2>&1"
    

# This rule creates a manifest file for the IsoQuant counts.
rule isoquant_counts:
    input:
        counts=expand(os.path.join(config["output_dir"], "isoquant","{{experiment}}","{condition}","{condition}.transcript_model_grouped_counts.tsv"), condition=grouped.keys()),
    output:
        manifest=os.path.join(config["output_dir"], "isoquant", "tama","{experiment}","count_manifest.tsv")
    conda:
        ENVS + "/kallisto_counts.yaml"
    run:
        with open(output.manifest, "w") as f:
            for condition, count_file in zip(grouped.keys(), input.counts):
                f.write(f"{condition}\t{count_file}\n")
