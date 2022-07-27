import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####


#WRAPPER_PREFIX='https://github.com/hans-vg/snakemake-wrappers/raw'
WRAPPER_PREFIX='https://raw.githubusercontent.com/hans-vg/snakemake-wrappers'

configfile: "config.yaml"

wildcard_constraints:
    sample="[\w-]+",
    trimmer="\w+"

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
#aligners=config["params"]["aligners"].split(",")
trimmers=config["params"]["trimmers"].split(",")
#print("Aligners:", aligners)
print("Trimmers:", trimmers)

cwd = os.getcwd()
print("Cwd:", cwd)

def strip_suffix(pattern, suffix):
    return pattern[: -len(suffix)]

wrappers_version="v1.7.1"

##### target rules #####
rule all:
    input:
        "qc/multiqc_report_pretrim.html",
        expand("trimmed/{trimmer}/{unit.sample}.{unit.unit}.1.fastq.gz", trimmer=trimmers, unit=units.itertuples()),
        expand("trimmed/{trimmer}/{unit.sample}.{unit.unit}.2.fastq.gz", trimmer=trimmers, unit=units.itertuples()),
        #expand("salmon/{trimmer}/{unit.sample}.{unit.unit}/lib_format_counts.json", trimmer=trimmers, unit=units.itertuples()),
        expand("salmon/{trimmer}/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}.lib_format_counts.json", trimmer=trimmers, unit=units.itertuples()),
        expand("qc/multiqc_report_salmon_{trimmer}.html", trimmer=trimmers),



include: "rules/qc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
