# This rule creates a manifest file for flair quantification.
rule flair_manifest:
    input:
        metadata=config["metadata"]
    output:
        manifest=os.path.join(config["output_dir"], "flair", "manifest","{experiment}","{condition}_manifest.tsv")
    conda:
        ENVS + "/kallisto_counts.yaml"
    run:
        # Create the new dataframe with the required columns
        output_df = metadata[['sample', 'fastq', 'condition']].copy()
        # filter the df to include only the condition of interest
        output_df = output_df[output_df['condition'] == wildcards.condition]
        output_df['batch'] = "batch1"

        # Reorder columns
        output_df = output_df[['sample', 'condition', 'batch', 'fastq']]

        # Save to TSV file
        output_tsv_path = output[0]
        output_df.to_csv(output_tsv_path, sep='\t', index=False, header=False)

rule merge_bam_by_condition:
    input:
        bam=lambda w: grouped[w.condition]["aligned"]
    output:
        merged_bam=temp(os.path.join(config["output_dir"], "flair", "{experiment}","{condition}_merged.bam")),
        merged_bam_index=temp(os.path.join(config["output_dir"], "flair", "{experiment}","{condition}_merged.bam.csi"))
    conda:
        ENVS + "/align.yaml"
    threads: 4
    resources:
        mem_mb=5000,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/flair/{experiment}/{condition}_merge_bam.log"
    benchmark:
        BENCHMARKS + "/flair/{experiment}/{condition}_merge_bam.txt"
    params:
        index="--write-index"
    shell:
        "samtools merge {params.index} -@ {threads} {output.merged_bam} {input.bam} > {log} 2>&1"

# use the new trnascriptome module in flair3 to generate the transcriptome
rule flair_transcriptome:
    input:
        bam=rules.merge_bam_by_condition.output.merged_bam,
        bam_index=rules.merge_bam_by_condition.output.merged_bam_index,
        genome=config["reference_genome"],
        annotation=config["reference_annotation"],
    output:
        gtf=os.path.join(config["output_dir"], "flair", "01_transcriptome", "{experiment}","{condition}.isoforms.gtf"),
        fasta=os.path.join(config["output_dir"], "flair", "01_transcriptome", "{experiment}","{condition}.isoforms.fa"),
        bed=os.path.join(config["output_dir"], "flair", "01_transcriptome", "{experiment}","{condition}.isoforms.bed")
    conda:
        ENVS + "/flair3.yaml"
    threads: 8
    resources:
        mem_mb=flair_transcriptome_memory,
        slurm_extra=flair_queue,
        runtime=flair_time
    log:
        LOGS + "/flair/{experiment}/{condition}_transcriptome.log"
    benchmark:
        BENCHMARKS + "/flair/{experiment}/{condition}_transcriptome.txt"
    params:
        prefix=os.path.join(config["output_dir"], "flair", "01_transcriptome", "{experiment}", "{condition}"),
        settings="--stringent --check_splice --generate_map --annotation_reliant generate"
    shell:
        """
        flair transcriptome --threads {threads} -b {input.bam} -g {input.genome} \
        --gtf {input.annotation} --output {params.prefix} {params.settings} \
         > {log} 2>&1
        """

# I am struggling with duplicated transcript_ids for the SIRVs, which can happen 
# in flair transcriptome. This rule renames those transcript_ids in the GTF file 
# to ensure they are unique.
rule deduplicate_gtf_names:
    input:
        gtf=rules.flair_transcriptome.output.gtf,
        script=SCRIPTS + "/dedup_gtf_names.py"
    output:
        gtf=os.path.join(config["output_dir"], "flair", "01_transcriptome", "{experiment}","{condition}.isoforms.dedup.gtf")
    conda:
        ENVS + "/flair3.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/flair/{experiment}/{condition}_dedup_gtf.log"
    benchmark:
        BENCHMARKS + "/flair/{experiment}/{condition}_dedup_gtf.txt"
    shell:
        """
        python {input.script} --input-gtf {input.gtf} --output-gtf {output.gtf} > {log} 2>&1
        """

# Generate transcript sequences from the deduplicated GTF and reference genome using gffread.
rule flair_generate_fasta:
    input:
        gtf=rules.deduplicate_gtf_names.output.gtf,
        genome=config["reference_genome"]
    output:
        fasta=os.path.join(config["output_dir"], "flair", "01_transcriptome", "{experiment}","{condition}.isoforms.dedup.fa")
    conda:
        ENVS + "/gffread.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/flair/{experiment}/{condition}_fasta_gen.log"
    benchmark:
        BENCHMARKS + "/flair/{experiment}/{condition}_fasta_gen.txt"
    shell:
        "gffread {input.gtf} -g {input.genome} -w {output.fasta} > {log} 2>&1"

# This rule quantifies the isoforms using flair quantify.
rule flair_quantify:
    input:
        manifest=rules.flair_manifest.output.manifest,
        fasta=rules.flair_generate_fasta.output.fasta,
    output:
        output=os.path.join(config["output_dir"], "flair", "02_quantify", "{experiment}","{condition}.counts.tsv")
    conda:
        ENVS + "/flair3.yaml"
    threads: 8
    resources:
        mem_mb=flair_memory,
        slurm_extra="'--qos=medium'",
        runtime="72h"
    log:
        LOGS + "/flair/{experiment}/{condition}_quantify.log"
    benchmark:
        BENCHMARKS + "/flair/{experiment}/{condition}_quantify.txt"
    params:
        prefix=os.path.join(config["output_dir"], "flair", "02_quantify", "{experiment}", "{condition}"),
        extra="--sample_id_only --quality 0 --generate_map",
        tmp_dir=os.path.join(config["output_dir"], "flair", "02_quantify", "{experiment}", "{condition}")
    shell:
        """
        flair quantify --threads {threads} --reads_manifest {input.manifest} \
            --isoforms {input.fasta} --output {params.prefix} \
            {params.extra} --temp_dir {params.tmp_dir} > {log} 2>&1
        """

# This rule creates a manifest file for the flair counts.
rule flair_counts:
    input:
        counts=expand(os.path.join(config["output_dir"], "flair","02_quantify","{{experiment}}","{condition}.counts.tsv"), condition=grouped.keys()),
    output:
        manifest=os.path.join(config["output_dir"], "flair", "tama","{experiment}","count_manifest.tsv")
    conda:
        ENVS + "/kallisto_counts.yaml"
    run:
        with open(output.manifest, "w") as f:
            for condition, count_file in zip(grouped.keys(), input.counts):
                f.write(f"{condition}\t{count_file}\n")

# This rule converts the flair GTF output to BED format.
rule flair_bed:
    input:
        gtf=rules.deduplicate_gtf_names.output.gtf,
        gtf2bed=rules.prepare_tama.output.gtf2bed,
        script=SCRIPTS + "/gtf2bed.sh"
    output:
        bed=os.path.join(config["output_dir"], "flair","tama","{experiment}","{condition}.bed")
    conda:
        ENVS + "/tama.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/flair/{experiment}/{condition}_flair_bed.log"
    benchmark:
        BENCHMARKS + "/flair/{experiment}/{condition}_flair_bed.txt"
    shell:
        "bash {input.script} {input.gtf} {output} {input.gtf2bed} > {log} 2>&1"