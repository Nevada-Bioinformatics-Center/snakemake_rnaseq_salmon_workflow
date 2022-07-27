
def get_trim_fastq1(wildcards):
    fq1 = expand("trimmed/{trimmer}/{sample}.{unit}.1.fastq.gz", **wildcards)
    return fq1

def get_trim_fastq2(wildcards):
    fq2 = expand("trimmed/{trimmer}/{sample}.{unit}.2.fastq.gz", **wildcards)
    return fq2

rule salmon_index:
    input:
        sequences=config["ref"]["transcriptomefa"],
    output:
        directory("salmon/transcriptome_index/"),
    log:
        "logs/salmon/transcriptome_index.log",
    threads: 16
    resources: time_min=320, mem_mb=20000, cpus=16
    params:
        # optional parameters
        extra="",
    wrapper:
        f"{wrappers_version}/bio/salmon/index"
        #"v1.7.1/bio/salmon/index"


rule salmon_quant_reads:
    input:
        r1="trimmed/{trimmer}/{sample}.{unit}.1.fastq.gz",
        r2="trimmed/{trimmer}/{sample}.{unit}.2.fastq.gz",
        index="salmon/transcriptome_index/",
    output:
        quant="salmon/{trimmer}/{sample}.{unit}/quant.sf",
        lib="salmon/{trimmer}/{sample}.{unit}/lib_format_counts.json",
    log:
        "logs/salmon/{trimmer}/{sample}.{unit}.log",
    params:
        # optional parameters
        libtype="A",
        extra="",
    threads: 16
    resources: time_min=320, mem_mb=20000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/salmon/quant"
        #"v1.7.1/bio/salmon/quant"

rule symlink_salmon_quants:
    input:
        quant="salmon/{trimmer}/{sample}.{unit}/quant.sf",
        lib="salmon/{trimmer}/{sample}.{unit}/lib_format_counts.json",
    output:
        quant="salmon/{trimmer}/{sample}.{unit}/{sample}.{unit}.quant.sf",
        lib="salmon/{trimmer}/{sample}.{unit}/{sample}.{unit}.lib_format_counts.json",
    params:
        #quant="{sample}.{unit}.quant.sf",
        #lib="{sample}.{unit}.lib_format_counts.json",
        quant="../../../salmon/{trimmer}/{sample}.{unit}/quant.sf",
        lib="../../../salmon/{trimmer}/{sample}.{unit}/lib_format_counts.json",
    threads: 1
    resources: time_min=320, mem_mb=2000, cpus=1
    shell:
        "ln -s {params.quant} {output.quant} && ln -s {params.lib} {output.lib}"
