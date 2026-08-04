# This rule creates a manifest file for TAMA merge.
rule tama_manifest:
    input:
        beds=expand(os.path.join(config["output_dir"], "{{pipeline}}","tama","{{experiment}}","{condition}.bed"), condition=grouped.keys()),
    output:
        manifest=os.path.join(config["output_dir"], "{pipeline}", "tama","{experiment}","tama_manifest.txt")
    conda:
        ENVS + "/kallisto_counts.yaml"
    run:
        with open(output.manifest, "w") as f:
            for condition, bed in zip(grouped.keys(), input.beds):
                f.write(f"{bed}\tcapped\t1,1,1\t{condition}\n")

# This rule merges transcript models from different conditions using TAMA.
rule tama_merge:
    input:
        manifest=rules.tama_manifest.output.manifest,
        tama_merge=rules.prepare_tama.output.tama_merge
    output:
        merged_bed=os.path.join(config["output_dir"], "{pipeline}","tama", "{experiment}","merged.bed"),
        gene_report=os.path.join(config["output_dir"], "{pipeline}","tama", "{experiment}","merged_gene_report.txt"),
        merge=os.path.join(config["output_dir"], "{pipeline}","tama", "{experiment}","merged_merge.txt"),
        trans_report=os.path.join(config["output_dir"], "{pipeline}","tama", "{experiment}","merged_trans_report.txt")
    conda:
        ENVS + "/tama.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/tama/{experiment}/{pipeline}/tama_merge.log"
    benchmark:
        BENCHMARKS + "/tama/{experiment}/{pipeline}/tama_merge.txt"
    params:
        prefix=os.path.join(config["output_dir"], "{pipeline}","tama", "{experiment}", "merged"),
        options="-m 0 -d merge_dup -a 50 -z 50"
    shell:
        "python {input.tama_merge} -f {input.manifest} -p {params} > {log} 2>&1"

# This rule creates a FASTA file from the merged BED file.
rule tama_merge_fasta:
    input:
        bed=rules.tama_merge.output.merged_bed,
        reference_fasta=config["reference_genome"]
    output:
        fasta=os.path.join(config["output_dir"], "{pipeline}","tama", "{experiment}","merged.fasta")
    conda:
        ENVS + "/gffread.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/tama/{experiment}/{pipeline}/tama_merge_fasta.log"
    benchmark:
        BENCHMARKS + "/tama/{experiment}/{pipeline}/tama_merge_fasta.txt"
    shell:
        "gffread -w {output} -g {input.reference_fasta} {input.bed} > {log} 2>&1"

# This rule calculates statistics for the merged transcriptome.
rule tama_transcriptome_stats:
    input:
        fasta=rules.tama_merge_fasta.output.fasta
    output:
        transcriptome_stats=os.path.join(config["output_dir"], "{pipeline}","tama", "{experiment}","{experiment}_transcript_models.tsv")
    conda:
        ENVS + "/seqkit.yaml"
    threads: 8
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/tama/{experiment}/{pipeline}/tama_merge_transcriptome_stats.log"
    benchmark:
        BENCHMARKS + "/tama/{experiment}/{pipeline}/tama_merge_transcriptome_stats.txt"
    shell:
        "(seqkit fx2tab -i -n -l -g -H {input.fasta} | cut -f2 -d';') > {output} 2> {log}"

# This rule merges the quantification results from different conditions.
rule tama_merge_counts:
    input:
        counts_manifest=os.path.join(config["output_dir"], "{pipeline}","tama","{experiment}","count_manifest.tsv"),
        tama_file=rules.tama_merge.output.merge,
        merge_script=SCRIPTS + "/merge_quantification.R"
    output:
        counts=os.path.join(config["output_dir"], "{pipeline}","tama","{experiment}","{experiment}_counts.tsv")
    conda:
        ENVS + "/tidy.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/tama/{experiment}/{pipeline}/tama_merge_counts.log"
    benchmark:
        BENCHMARKS + "/tama/{experiment}/{pipeline}/tama_merge_counts.txt"
    shell:
        "Rscript {input.merge_script} {input.counts_manifest} {input.tama_file} {output} > {log} 2>&1"


# This rule filters the merged BED file based on the quantified transcripts.
rule tama_filter_bed:
    input:
        bed=rules.tama_merge.output.merged_bed,
        counts=rules.tama_merge_counts.output.counts,
    output:
        filtered_bed=os.path.join(config["output_dir"], "{pipeline}","tama", "{experiment}","filtered.bed")
    conda:
        ENVS + "/kallisto_counts.yaml" # Anything with pandas should do the trick
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/tama/{experiment}/{pipeline}/tama_filter_bed.log"
    benchmark:
        BENCHMARKS + "/tama/{experiment}/{pipeline}/tama_filter_bed.txt"
    run:
        import pandas as pd
        df = pd.read_csv(input.counts, sep="\t", header=0, index_col=0)
        trasncripts_id = set(df.index)
        with open(input.bed, "r") as f:
            with open(output.filtered_bed, "w") as out:
                for line in f:
                    if line.startswith("#"):
                        out.write(line)
                    else:
                        transcript_id = line.split("\t")[3].split(";")[1] # tama transcript ids are formated like gene_id;transcript_id 
                        if transcript_id in trasncripts_id:
                            out.write(line)

# This rule converts the filtered BED file to GTF format.
rule tama_merge_gtf:
    input:
        bed=rules.tama_filter_bed.output.filtered_bed,
        bed2gtf=rules.prepare_tama.output.bed2gtf,
    output:
        gtf=os.path.join(config["output_dir"], "{pipeline}","tama", "{experiment}","merged.gtf")
    conda:
        ENVS + "/tama.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/tama/{experiment}/{pipeline}/tama_merge_gtf.log"
    benchmark:
        BENCHMARKS + "/tama/{experiment}/{pipeline}/tama_merge_gtf.txt"
    shell:
        "python {input.bed2gtf} {input.bed} {output.gtf} > {log} 2>&1"


# This rule performs SQANTI3 QC on the merged GTF file.
rule sqanti3_merge:
    input:
        gtf=rules.tama_merge_gtf.output.gtf,
        sqanti3_script=rules.prepare_sqanti.output.sqanti_qc,
        reference_annotation=config["reference_annotation"],
        reference_genome=config["reference_genome"],
    output:
        sqanti3_classification=os.path.join(config["output_dir"], "{pipeline}","tama", "{experiment}","sqanti3", "merged_classification.txt"),
        out_dir=directory(os.path.join(config["output_dir"], "{pipeline}","tama", "{experiment}","sqanti3")),
    conda:
        ENVS + "/SQANTI3.yml"
    threads: 2
    resources:
        mem_mb=sqanti_memory,
        slurm_extra=sqanti_queue,
        runtime=sqanti_time
    log:
        LOGS + "/tama/{experiment}/{pipeline}/sqanti3_merge.log"
    benchmark:
        BENCHMARKS + "/tama/{experiment}/{pipeline}/sqanti3_merge.txt"
    shell:
        '''
        python3 {input.sqanti3_script} {input.gtf} {input.reference_annotation} \
            {input.reference_genome} --dir {output.out_dir} --report skip \
            --skipORF --output merged -t {threads} > {log} 2>&1
        '''    

# This rule assigns reference transcripts to the merged transcripts.
rule assign_reference:
    input:
        sqanti3_classification=rules.sqanti3_merge.output.sqanti3_classification,
        quantification=rules.tama_merge_counts.output.counts,
        script=SCRIPTS + "/assign_reference.py",
        transcriptome_stats=rules.tama_transcriptome_stats.output.transcriptome_stats,
    output:
        condensed_class=os.path.join(config["output_dir"], "{pipeline}","NOIseq", "{experiment}_classification.txt"),
        condensed_counts=os.path.join(config["output_dir"], "{pipeline}","NOIseq", "{experiment}_counts.tsv"),
        trans_map=os.path.join(config["output_dir"], "{pipeline}","NOIseq", "{experiment}_collapse_map.csv"),
        condensed_categories=os.path.join(config["output_dir"], "{pipeline}","NOIseq", "{experiment}_structural_category.tsv"),
        condensed_stats=os.path.join(config["output_dir"], "{pipeline}","NOIseq", "{experiment}_transcript_models.tsv"),
        condensed_transcript_to_gene=os.path.join(config["output_dir"], "{pipeline}","NOIseq", "{experiment}_transcript_to_gene.tsv"),
    conda:
        ENVS + "/kallisto_counts.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/tama/{experiment}/{pipeline}/assign_reference.log"
    benchmark:
        BENCHMARKS + "/tama/{experiment}/{pipeline}/assign_reference.txt"
    params:
        out_dir=os.path.join(config["output_dir"], "{pipeline}","NOIseq", "{experiment}"),
    shell:
        '''
        python3 {input.script} -c {input.sqanti3_classification} -q {input.quantification} \
            -m {input.transcriptome_stats} -o {params.out_dir} > {log} 2>&1
        cut -f1,6 {output.condensed_class} > {output.condensed_categories} 2>> {log}
        '''