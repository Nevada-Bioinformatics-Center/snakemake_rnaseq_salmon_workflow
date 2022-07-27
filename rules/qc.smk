def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

def get_fastq1(wildcards):
    fq1 = units.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna()
    return fq1

def get_fastq2(wildcards):
    fq2 = units.loc[(wildcards.sample, wildcards.unit), ["fq2"]].dropna()
    return fq2


rule fastqc_pretrim_r1:
    input:
       get_fastq1
    output:
        html="qc/fastqc_pretrim/{trimmer}/{sample}.{unit}_r1.html",
        zip="qc/fastqc_pretrim/{trimmer}/{sample}.{unit}_r1_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_pretrim/{trimmer}/{sample}.{unit}_r1.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    threads: 1
    wrapper:
        #"v0.75.0/bio/fastqc"
        f"{wrappers_version}/bio/fastqc"

rule fastqc_pretrim_r2:
    input:
       get_fastq2
    output:
        html="qc/fastqc_pretrim/{trimmer}/{sample}.{unit}_r2.html",
        zip="qc/fastqc_pretrim/{trimmer}/{sample}.{unit}_r2_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_pretrim/{trimmer}/{sample}.{unit}_r2.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    threads: 1
    wrapper:
        #"v0.75.0/bio/fastqc"
        f"{wrappers_version}/bio/fastqc"

rule fastqc_posttrim_r1:
    input:
        "trimmed/{trimmer}/{sample}.{unit}.1.fastq.gz"
    output:
        html="qc/fastqc_posttrim/{trimmer}/{sample}.{unit}_r1.html",
        zip="qc/fastqc_posttrim/{trimmer}/{sample}.{unit}_r1_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_posttrim/{trimmer}/{sample}.{unit}_r1.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    threads: 1
    wrapper:
        f"{wrappers_version}/bio/fastqc"

rule fastqc_posttrim_r2:
    input:
        "trimmed/{trimmer}/{sample}.{unit}.2.fastq.gz"
    output:
        html="qc/fastqc_posttrim/{trimmer}/{sample}.{unit}_r2.html",
        zip="qc/fastqc_posttrim/{trimmer}/{sample}.{unit}_r2_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_posttrim/{trimmer}/{sample}.{unit}_r2.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    threads: 1
    wrapper:
        #"v0.75.0/bio/fastqc"
        f"{wrappers_version}/bio/fastqc"

rule multiqc_pre:
    input:
        expand("qc/fastqc_pretrim/{trimmer}/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples(), trimmer=trimmers),
        expand("qc/fastqc_pretrim/{trimmer}/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples(), trimmer=trimmers)
    output:
        "qc/multiqc_report_pretrim.html"
    log:
        "logs/multiqc_pre.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        #"0.84.0/bio/multiqc"
        f"{wrappers_version}/bio/multiqc"
        #f"{WRAPPER_PREFIX}/master/bio/multiqc"

rule multiqc_post_trimmomatic:
    input:
        expand("logs/trimmomatic/{unit.sample}.{unit.unit}.log", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimmomatic/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimmomatic/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_posttrim_trimmomatic.html"
    log:
        "logs/multiqc_posttrim_trimmomatic.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        #"0.84.0/bio/multiqc"
        f"{wrappers_version}/bio/multiqc"
        #f"{WRAPPER_PREFIX}/master/bio/multiqc"

rule multiqc_post_fastp:
    input:
        expand("report/fastp/{unit.sample}.{unit.unit}.fastp.json", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_posttrim_fastp.html"
    log:
        "logs/multiqc_posttrim_fastp.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

rule multiqc_post_trimgalore:
    input:
        expand("logs/trimgalore/{unit.sample}.{unit.unit}.log", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimgalore/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimgalore/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_posttrim_trimgalore.html"
    log:
        "logs/multiqc_posttrim_trimgalore.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        #"0.84.0/bio/multiqc"
        f"{wrappers_version}/bio/multiqc"
        #f"{WRAPPER_PREFIX}/master/bio/multiqc"

rule multiqc_salmon_trimmomatic:
    input:
        expand("salmon/trimmomatic/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}.lib_format_counts.json", unit=units.itertuples()),
        expand("logs/trimmomatic/{unit.sample}.{unit.unit}.log", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimmomatic/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimmomatic/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_salmon_trimmomatic.html"
    log:
        "logs/multiqc_salmon_trimmomatic.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

rule multiqc_salmon_fastp:
    input:
        expand("salmon/fastp/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}.lib_format_counts.json", unit=units.itertuples()),
        expand("report/fastp/{unit.sample}.{unit.unit}.fastp.json", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_salmon_fastp.html"
    log:
        "logs/multiqc_salmon_fastp.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"


