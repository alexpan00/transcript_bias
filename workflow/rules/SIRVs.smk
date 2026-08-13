# This rule counts the number of reads that map to the SIRV transcripts.
rule count_SIRV_reads:
    input:
        bam=lambda wildcards: grouped_by_sample[wildcards.sample]["aligned"]
    output:
        counts=os.path.join(config["output_dir"], "SIRVs", "{sample}", "SIRV_counts.tsv")
    conda:
        ENVS + "/align.yaml"
    threads: 4
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/SIRVs/counts/{sample}.log"
    benchmark:
        BENCHMARKS + "/SIRVs/counts/{sample}.txt"
    shell:
        '''
        samtools idxstats {input.bam} | cut -f 1,3 | grep SIRV | awk -F "\t" '{{sum += $2}} END {{print sum}}' > {output.counts} 2> {log}
        '''

# This rule calculates the expected SIRV counts based on the total number of reads.
rule SIRV_expected_quant:
    input:
        counts=rules.count_SIRV_reads.output.counts,
        script=SCRIPTS + "/quantification/SIRV_expected_quant.py",
        mix_matrix=config.get("sirv_mixes", "config/sirv_mixes.json")
    output:
        expected=os.path.join(config["output_dir"], "SIRVs", "{sample}", "SIRV_expected.tsv")
    conda:
        ENVS + "/kallisto_counts.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/SIRVs/expected/{sample}.log"
    benchmark:
        BENCHMARKS + "/SIRVs/expected/{sample}.txt"
    params:
        sirv_set=lambda wildcards: grouped_by_sample[wildcards.sample]["SIRV"],
    shell:
        '''
        python {input.script} --reads {input.counts} --set {params.sirv_set} \
            --mix_matrix {input.mix_matrix} --out {output.expected} > {log} 2>&1
        '''


# This rule merges the expected SIRV counts from all samples into a single file.
rule merge_SIRV_counts:
    input:
        counts=expand(os.path.join(config["output_dir"], "SIRVs", "{sample}", "SIRV_expected.tsv"), sample=list(grouped_by_sample.keys()))
    output:
        counts=os.path.join(config["output_dir"], "SIRVs", "SIRV_counts.tsv")
    conda:
        ENVS + "/kallisto_counts.yaml"
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/SIRVs/merge_counts.log"
    benchmark:
        BENCHMARKS + "/SIRVs/merge_counts.txt"
    run:
        import pandas as pd
        import os
        try:
            df = pd.read_csv(input.counts[0], sep="\t", index_col=0)
            # replace the exp_count column by the sample name
            sample_name = os.path.basename(os.path.dirname(input.counts[0]))
            df.rename(columns={'exp_count': sample_name}, inplace=True)
            print(df)
            # loop over the other files and create a new column for each with the sample name
            for input_file in input.counts[1:]:
                sample_name = os.path.basename(os.path.dirname(input_file))
                df2 = pd.read_csv(input_file, sep='\t', index_col=0)
                print(df2)
                df[sample_name] = df2['exp_count']
            
            # write the output file
            df.to_csv(output.counts, sep='\t', index=True, index_label='SIRV')
        except Exception as e:
            with open(log[0], "w") as f:
                f.write(f"Error: {e}\n")

# This rule prepares the SIRV counts for NOISeq analysis.
rule SIRV_counts_NOISeq:
    input:
        counts=rules.merge_SIRV_counts.output.counts,
        script=SCRIPTS + "/normalization/SIRV_NOIseq.R",
        factors=config["factors"],
        sirvs_info=config.get("sirvs_info", "sirvs_info.csv")
    output:
        noiseqcounts=os.path.join(config["output_dir"], "SIRVs", "NOIseq", "{experiment}_SIRV_counts.rds")
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/SIRVs/NOIseq/{experiment}.log"
    benchmark:
        BENCHMARKS + "/SIRVs/NOIseq/{experiment}.txt"
    conda:
        ENVS + "/NOIseq.yaml"
    shell:
        '''
        Rscript {input.script} {input.counts} {input.factors} {input.sirvs_info} \
            {output.noiseqcounts} > {log} 2>&1
        '''


# This rule prepares the ERCC counts for NOISeq analysis.
rule ERCC_counts_NOISeq:
    input:
        counts=config.get("ERCC_counts", ""),
        script=SCRIPTS + "/normalization/ERCC_NOIseq.R",
        factors=config["factors"]
    output:
        noiseqcounts=os.path.join(config["output_dir"], "SIRVs", "NOIseq", "{experiment}_ERCC_counts.rds")
    threads: 1
    resources:
        mem_mb=generic_memory,
        slurm_extra="'--qos=short'"
    log:
        LOGS + "/ERCC/NOIseq/{experiment}.log"
    benchmark:
        BENCHMARKS + "/ERCC/NOIseq/{experiment}.txt"
    conda:
        ENVS + "/NOIseq.yaml"
    shell:
        '''
        Rscript {input.script} {input.counts} {input.factors} \
            {output.noiseqcounts} > {log} 2>&1
        '''