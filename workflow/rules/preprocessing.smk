# This rule indexes the BAM files.
rule index_bam:
    input:
        bam="{bam_path}/{sample}.bam"
    output:
        index="{bam_path}/{sample}.bam.bai"
    conda:
        ENVS + "/align.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/index_bam/{bam_path}/{sample}_index.log"
    benchmark:
        BENCHMARKS + "/index_bam/{bam_path}/{sample}_index.txt"
    shell:
        "samtools index {input.bam} {output} > {log} 2>&1"

# This rule prepares the SQANTI3 environment by cloning the git repository.
rule prepare_sqanti:
    output:
        sqanti_qc=os.path.join(config["sqanti_dir"], "sqanti3_qc.py"),
        sqanti_filter=os.path.join(config["sqanti_dir"], "sqanti3_filter.py"),
        sqanti_rescue=os.path.join(config["sqanti_dir"], "sqanti3_rescue.py")
    conda:
        ENVS + "/SQANTI3.yml"
    threads: 1
    log:    
        LOGS + "/sq3/install.log"
    shell:
        "rm -rf {sqanti_dir} ; "
        "git clone https://github.com/alexpan00/SQANTI3.git {sqanti_dir} 2> {log}"

# This rule prepares the TAMA environment by cloning the git repository.
rule prepare_tama:
    output:
        tama_merge=os.path.join(tama_dir, "tama_merge.py"),
        gtf2bed=os.path.join(tama_dir, "tama_go/format_converter/tama_format_gff_to_bed12_cupcake.py"),
        bed2gtf=os.path.join(tama_dir, "tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py"),
    conda:
        ENVS + "/tama.yaml"
    threads: 1
    log:
        LOGS + "/tama/install.log"
    shell:
        "rm -rf {tama_dir} ; "
        "git clone https://github.com/GenomeRIK/tama.git {tama_dir} 2> {log}"

# Prepare transcript to gene mapping file based on the reference annotation
rule prepare_transcript_to_gene_map:
    input:
        annotation=config["reference_annotation"]
    output:
        transcript_to_gene=os.path.join(config["output_dir"], "transcript_to_gene.tsv")
    conda:
        ENVS + "/gffread.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/prepare_transcript_to_gene_map.log"
    shell:
        """
        gffread {input.annotation} --table @id,@geneid -o {output.transcript_to_gene} 2> {log}
        """
# This rule prepares a file that links samples to their corresponding files basenames
rule prepare_file_to_sample:
    input:
        metadata=config["metadata"]
    output:
        metadata_extended=os.path.join(os.path.join(config["output_dir"], "metadata_extended.tsv"))
    conda:
        ENVS + "/tama.yaml"
    threads: 1
    log:
        LOGS + "/prepare_file_to_sample.log"
    run:
        from os.path import getsize, basename

        def add_basenames(df):
            """Add bam_basename and fastq_basename columns to the metadata DataFrame."""
            df["bam_basename"] = df["aligned"].apply(lambda x: basename(x).split('.')[0])
            df["fastq_basename"] = df["fastq"].apply(lambda x: basename(x).split('.')[0])
            return df
        metadata_extended = add_basenames(metadata.copy())
        # select sample and basenames columns
        df = metadata_extended[["sample", "bam_basename", "fastq_basename"]]
        df.to_csv(output.metadata_extended, sep="\t", index=False)


# This rule copies the workflow configuration, metadata, and factors into a hidden .run_metadata folder in output_dir
rule save_run_metadata:
    input:
        metadata=config["metadata"],
        factors=config["factors"] if config.get("factors") else []
    output:
        meta_dir=directory(os.path.join(config["output_dir"], ".run_metadata")),
        metadata_copy=os.path.join(config["output_dir"], ".run_metadata", "metadata.csv"),
        config_copy=os.path.join(config["output_dir"], ".run_metadata", "config.yml")
    params:
        factors=config.get("factors", "")
    run:
        import shutil
        import yaml
        
        os.makedirs(output.meta_dir, exist_ok=True)
        shutil.copy(input.metadata, output.metadata_copy)
        
        if params.factors and os.path.exists(params.factors):
            shutil.copy(params.factors, os.path.join(output.meta_dir, "factors.csv"))
            
        with open(output.config_copy, "w") as f:
            yaml.dump(dict(config), f, default_flow_style=False)


