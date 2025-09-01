#!/usr/bin/env nextflow

params.reads1 = ""
params.reads2 = ""
params.sample = "sample"

process FASTP {
    tag "$params.sample"
    publishDir "results", mode: "copy"

    input:
    path reads1 from params.reads1
    path reads2 from params.reads2

    output:
    path "${params.sample}.trimmed_merged.fastq.gz"
    path "fastp_report.html"
    path "fastp_report.json"

    script:
    """
    fastp \
        --in1 $reads1 \
        --in2 $reads2 \
        -p -c --merge \
        --merged_out ${params.sample}.trimmed_merged.fastq.gz \
        -h fastp_report.html \
        -j fastp_report.json \
        -w 4 \
        -l 30
    """
}

