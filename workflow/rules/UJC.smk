# Rule to get the ids of the transcripts after the counts filtering
rule whitelist_transcripts:
    input:
        NOISeq_object=rules.NOIseq_object.output.noiseq_obj,
        script=os.path.join(SCRIPTS, "get_transcripts_ids.R")
    output:
        transcripts_ids=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{experiment}_transcripts_ids.txt")
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=NOISeq_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{tool}/{experiment}_transcripts_ids.log"
    benchmark:
        BENCHMARKS + "/{tool}/{experiment}_transcripts_ids.txt"
    shell:
        '''
        Rscript {input.script} {input.NOISeq_object} {output.transcripts_ids} > {log} 2>&1
        '''

# use gffread to filter the gtf file with the ids of the transcripts
# and convert to table format
rule preprocess_gtf_transcript_ids:
    input:
        gtf=os.path.join(config["output_dir"], "{tool}","tama", "{experiment}","merged.gtf"),
        collapse_map=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{experiment}_collapse_map.csv"),
        script=os.path.join(SCRIPTS, "remap_gtf_transcript_ids.py")
    output:
        preprocessed_gtf=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{experiment}_remapped.gtf")
    conda:
        ENVS + "/isoquant.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{tool}/{experiment}_preprocess_gtf.log"
    benchmark:
        BENCHMARKS + "/{tool}/{experiment}_preprocess_gtf.txt"
    shell:
        '''
        python {input.script} {input.gtf} {input.collapse_map} {output.preprocessed_gtf} > {log} 2>&1
        '''

# use pipeline specific gtf for transcript reconstruction methods and reference
# annotation for quantification methods. Then filter the gtf with the ids of the transcripts
# detected using each tool.
rule filter_gtf_transcripts:
    input:
        transcripts_ids=rules.whitelist_transcripts.output.transcripts_ids,
        gtf=lambda wildcards: config["reference_annotation"] if wildcards.tool not in pipelines else rules.preprocess_gtf_transcript_ids.output.preprocessed_gtf
    output:
        filtered_tlf=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{experiment}_filtered.tlf")
    conda:
        ENVS + "/gffread.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{tool}/{experiment}_filter_gtf.log"
    benchmark:
        BENCHMARKS + "/{tool}/{experiment}_filter_gtf.txt"
    shell:
        '''
        gffread --ids {input.transcripts_ids} --tlf {input.gtf} -o {output.filtered_tlf} > {log} 2>&1
        '''

# generate the UJC results from the filtered tlf file
rule generate_UJC_results:
    input:
        transcripts=rules.filter_gtf_transcripts.output.filtered_tlf,
        script=os.path.join(SCRIPTS, "get_ujc_from_tlf.py")
    output:
        ujc=os.path.join(config["output_dir"], "{tool}", "NOIseq", "{experiment}_UJC.tsv")
    conda:
        ENVS + "/isoquant.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{tool}/{experiment}_ujc.log"
    benchmark:
        BENCHMARKS + "/{tool}/{experiment}_ujc.txt"
    shell:
        '''
        python {input.script} {input.transcripts} {output.ujc} > {log} 2>&1
        '''


# collect all UJC files and summarize their overlap
rule generate_UJC_upset_plot:
    input:
        ujcs=expand(os.path.join(config["output_dir"], "{tool}", "NOIseq", "{{experiment}}_UJC.tsv"), tool=user_tools),
        script=os.path.join(SCRIPTS, "plot_ujc_upset.R")
    output:
        plot=os.path.join(config["output_dir"], "{experiment}_UJC_upset.png"),
        table=os.path.join(config["output_dir"], "{experiment}_UJC_summary.tsv")
    conda:
        ENVS + "/NOIseq.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/{experiment}_ujc_upset.log"
    benchmark:
        BENCHMARKS + "/{experiment}_ujc_upset.txt"
    shell:
        '''
        Rscript {input.script} {input.ujcs} {output.plot} {output.table} > {log} 2>&1
        '''