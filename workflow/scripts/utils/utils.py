import os
from os.path import getsize
from math import log

BYTE2MB = 1/1024**2
H2M = 1/60
MIN2DAY = 1/(60*24)


def generic_memory(wc, input, attempt):
    """Calculate memory requirements for a generic job."""
    return 1.25*input.size_mb + 5*(attempt+1)*1024

def flair_transcriptome_memory(wc, input, attempt):
    """Calculate memory requirements for a flair transcriptome job."""
    return 0.3*generic_memory(wc, input, attempt)

def extra_memory(wc, input, attempt):
    """Calculate memory requirements for a heavy job."""
    return 4*input.size_mb + 10*(attempt+1)*1024

def bambu_memory(wc, input, attempt):
    """Calculate memory requirements for a bambu job."""
    return 0.5*extra_memory(wc, input, attempt)

def NOISeq_memory(wc, input):
    """Calculate memory requirements for a NOISeq job."""
    return 1.25*input.size_mb + 10*1024

def align_memory(wc, input):
    """Calculate memory requirements for an alignment job."""
    target_file = input[1] if len(input) > 1 else input[0]
    return min(max(15.5*getsize(target_file)*BYTE2MB, 5*1024), 1000*1024)

def index_memory(wc, input):
    """Calculate memory requirements for an indexing job."""
    return 1.25*input.size_mb + 8*1024


def flair_memory(wc, input, attempt):
    """Calculate memory requirements for a flair job."""
    return 1.25*input.size_mb + 15*1024*(attempt + 1)

def flair_time(wc, attempt):
    """Calculate time requirements for a flair job."""
    return 48/H2M + 12*(attempt - 1)/H2M

def flair_queue(wc, input):
    """get flair qos."""
    return "'--qos=medium'"

def sqanti_memory(wc, input):
    """Calculate memory requirements for a SQANTI3 job."""
    return 5.5*input.size_mb + 5*1024

def sqanti_time(wc, input):
    """Calculate time requirements for a SQANTI3 job."""
    if not input or len(input) == 0:
        return 6 / H2M
    file_path = str(input[0])
    if not os.path.exists(file_path):
        return 6 / H2M
    size_mb = max(getsize(file_path) * BYTE2MB, 1e-6)
    return max(6 / H2M, (-1 + 400 * log(size_mb)))

def sqanti_queue(wc, input):
    """Determine the queue for a SQANTI3 job based on its estimated runtime."""
    return "'--qos=medium'" if sqanti_time(wc, input) > 24/H2M else "'--qos=short'"

def group_by_condition(df):
    """Group a pandas DataFrame by the 'condition' column."""
    result = df.groupby("condition").apply(lambda g: {col: g[col].tolist() for col in df.columns if col != "condition"}).to_dict()
    return result

def group_by_sample(df):
    """Group a pandas DataFrame by the 'sample' column."""
    result = df.groupby("sample").apply(lambda g: {col: g[col].tolist() for col in df.columns if col != "sample"}).to_dict()
    return result

def format_sjtab(wildcards, sample_info):
    """Format the SJ.tab file paths for the SQANTI3 command."""
    x = sample_info[wildcards.condition].get("SJtab", "")
    return "-c " + ",".join(x) if x else ""

def check_tools(user_tools, tools):
    """Check if the user-provided tools are valid."""
    for tool in user_tools:
        if tool not in tools:
            raise ValueError(f"Tool '{tool}' not implemented. Available tools are: {', '.join(tools)}")

def validate_samples_and_factors(metadata_path, factors_path, report_factors=None):
    """
    Validates that metadata and factors files exist, contain required columns,
    and that all metadata samples are present in the factors file.
    """
    import pandas as pd
    if not os.path.exists(metadata_path):
        raise FileNotFoundError(f"Metadata file not found at: {metadata_path}")
    if not os.path.exists(factors_path):
        raise FileNotFoundError(f"Factors file not found at: {factors_path}")

    metadata_df = pd.read_csv(metadata_path)
    factors_df = pd.read_csv(factors_path)

    if "sample" not in metadata_df.columns:
        raise ValueError(f"Required 'sample' column missing from metadata file: {metadata_path}")
    
    factor_sample_col = "sample" if "sample" in factors_df.columns else ("Sample" if "Sample" in factors_df.columns else factors_df.columns[0])

    meta_samples = set(metadata_df["sample"].dropna().astype(str).str.strip())
    factor_samples = set(factors_df[factor_sample_col].dropna().astype(str).str.strip())

    meta_dups = metadata_df["sample"][metadata_df["sample"].duplicated()].tolist()
    if meta_dups:
        raise ValueError(f"Duplicate sample IDs found in metadata file ({metadata_path}): {meta_dups}")

    missing_in_factors = meta_samples - factor_samples
    if missing_in_factors:
        raise ValueError(f"Sample(s) present in metadata ({metadata_path}) but missing from factors ({factors_path}): {sorted(list(missing_in_factors))}")

    # Validate report_factors if specified
    if report_factors:
        missing_factors = []
        for factor in report_factors:
            if factor not in factors_df.columns:
                parts = factor.split("_")
                if not (len(parts) > 1 and all(p in factors_df.columns for p in parts)):
                    missing_factors.append(factor)
        if missing_factors:
            raise ValueError(f"The following factor(s) specified in 'report_factors' were not found in factors file ({factors_path}): {missing_factors}\nAvailable columns in factors file: {list(factors_df.columns)}")

def get_rule_resource(config, rule_name, resource_key, default_val_or_func):
    """
    Returns a resource-resolver function usable directly inside Snakemake rule directives.
    If config['resources'][rule_name][resource_key] is specified, returns that custom value.
    Otherwise falls back to default_val_or_func (which can be a function or a value).
    For memory resources (e.g., 'mem_mb', 'mem_mib', 'memory'), increases memory by 10% on each retry attempt.
    """
    def _resource_resolver(wildcards=None, input=None, attempt=1, *args, **kwargs):
        resources_config = config.get("resources", {})
        val = None
        if isinstance(resources_config, dict):
            rule_res = resources_config.get(rule_name, {})
            if isinstance(rule_res, dict) and resource_key in rule_res:
                val = rule_res[resource_key]
        if val is None:
            if callable(default_val_or_func):
                try:
                    val = default_val_or_func(wildcards, input, attempt)
                except TypeError:
                    try:
                        val = default_val_or_func(wildcards, input)
                    except TypeError:
                        try:
                            val = default_val_or_func(wildcards)
                        except TypeError:
                            val = default_val_or_func()
            else:
                val = default_val_or_func

        if (resource_key in ("mem_mb", "mem_mib", "memory") or resource_key.startswith("mem")) and isinstance(val, (int, float)):
            multiplier = 1.0 + 0.1 * (attempt - 1)
            val = int(val * multiplier)

        return val

    return _resource_resolver