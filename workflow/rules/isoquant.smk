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
        LOGS + "/isoquant_db/{experiment}_isoquant_db.log"
    benchmark:
        BENCHMARKS + "/isoquant_db/{experiment}_isoquant_db.txt"
    script:
        SCRIPTS + "/preprocessing/isoquant_db.py"


# This rule runs IsoQuant to perform transcript reconstruction and quantification.
rule isoquant:
    input:
        bam=lambda wildcards: grouped[wildcards.condition]["aligned"],
        index=lambda wildcards: [bam + ".bai" for bam in grouped[wildcards.condition]["aligned"]],
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
        prefix="--prefix {condition}"
    shell:
        "isoquant.py -t {threads} --bam {input.bam} {params} > {log} 2>&1"

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
