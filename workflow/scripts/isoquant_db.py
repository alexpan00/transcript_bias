import sys
from snakemake.script import snakemake
import gffutils

try:
    gffutils.create_db(
        snakemake.input.annotation,
        snakemake.output.db,
        force=True,
        keep_order=True,
        merge_strategy='error',
        sort_attribute_values=True,
        disable_infer_transcripts=False,
        disable_infer_genes=False
    )
except Exception as e:
    log_path = snakemake.log[0] if isinstance(snakemake.log, list) else str(snakemake.log)
    if log_path:
        with open(log_path, "w") as f:
            f.write(str(e))
    sys.exit(1)