# This rule runs bambu to perform transcript reconstruction and quantification.
rule bambu:
    input:
        bam=lambda wildcards: grouped[wildcards.condition]["aligned"],
        script=SCRIPTS + "/quantification/bambu.R"
    output:
        bambu_gtf=os.path.join(config["output_dir"], "bambu","{experiment}","{condition}.gtf"),
        bambu_counts=os.path.join(config["output_dir"], "bambu","{experiment}","{condition}_counts.tsv")
    conda:
        ENVS + "/bambu.yaml"
    threads: get_rule_resource(config, "bambu", "threads", 8)
    resources:
        mem_mb=get_rule_resource(config, "bambu", "mem_mb", bambu_memory),
        slurm_extra=get_rule_resource(config, "bambu", "slurm_extra", "'--qos=short'")
    log:
        LOGS + "/bambu/{experiment}/{condition}_bambu.log"
    benchmark:
        BENCHMARKS + "/bambu/{experiment}/{condition}_bambu.txt"
    params:
        genome=f"{config['reference_genome']}",
        annotation=f"{config['reference_annotation']}",
        output_dir=os.path.join(config["output_dir"], "bambu","{experiment}","{condition}"),
        bam_list=lambda wildcards: ",".join(grouped[wildcards.condition]["aligned"]),
        seed=config.get("seed", 42)
    shell:
        "Rscript {input.script} {params.genome} {params.annotation} {params.bam_list} {params.output_dir} {threads} {params.seed} > {log} 2>&1"

# this rule fixes bambu GTF formatting issues for downstream processing
rule fix_bambu_gtf:
    input:
        gtf=rules.bambu.output.bambu_gtf,
        script=SCRIPTS + "/preprocessing/fix_bambu_gtf.py"
    output:
        fixed_gtf=os.path.join(config["output_dir"], "bambu","{experiment}","{condition}_fixed.gtf")
    conda:
        ENVS + "/bambu.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/bambu/{experiment}/{condition}_fix_bambu_gtf.log"
    benchmark:
        BENCHMARKS + "/bambu/{experiment}/{condition}_fix_bambu_gtf.txt"
    shell:
        """
        python {input.script} {input.gtf} {output.fixed_gtf} > {log} 2>&1
        """


# This rule converts the bambu GTF output to BED format.
rule bambu_bed:
    input:
        gtf=rules.fix_bambu_gtf.output.fixed_gtf,
        gtf2bed=rules.prepare_tama.output.gtf2bed,
        script=SCRIPTS + "/preprocessing/gtf2bed.sh",
    output:
        bed=os.path.join(config["output_dir"], "bambu","tama","{experiment}","{condition}.bed")
    conda:
        ENVS + "/tama.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/bambu/{experiment}/{condition}_bambu_bed.log"
    benchmark:
        BENCHMARKS + "/bambu/{experiment}/{condition}_bambu_bed.txt"
    shell:
        """
        bash {input.script} {input.gtf} {output} {input.gtf2bed} > {log} 2>&1
        """

# This rule creates a manifest file for the bambu counts.
rule bambu_counts:
    input:
        counts=expand(os.path.join(config["output_dir"], "bambu","{{experiment}}","{condition}_counts.tsv"), condition=grouped.keys()),
    output:
        manifest=os.path.join(config["output_dir"], "bambu", "tama","{experiment}","count_manifest.tsv")
    conda:
        ENVS + "/kallisto_counts.yaml"
    run:
        with open(output.manifest, "w") as f:
            for condition, count_file in zip(grouped.keys(), input.counts):
                f.write(f"{condition}\t{count_file}\n")
