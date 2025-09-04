#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.samplesheet = ""

process FASTP {
    tag "$sample"
    publishDir "results", mode: "copy"

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path "${sample}.trimmed_merged.fastq.gz"
    path "${sample}.fastp_report.html"
    path "${sample}.fastp_report.json"

    script:
    """
    fastp \
        --in1 $read1 \
        --in2 $read2 \
        -p -c --merge \
        --merged_out ${sample}.trimmed_merged.fastq.gz \
        -h ${sample}.fastp_report.html \
        -j ${sample}.fastp_report.json \
        -w 2 \
        -l 30
    """
}

workflow {
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true, sep:'\t')
        .map { row ->
            tuple(row.sample, file(row.read1), file(row.read2))
        }
        | FASTP
}

