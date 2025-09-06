#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.samplesheet = ""

params.db_mammalian   = "DBDIR_VERTEBRATE_MAMMALIAN"
params.db_other       = "DBDIR_VERTEBRATE_OTHER"
params.db_plant       = "DBDIR_PLANT"
params.db_invertebrate= "DBDIR_INVERTEBRATE"
params.db_microbe     = "DBDIR_MICROBE"

process FASTP {
    tag "$sample"
    publishDir "results/${sample}/fastp", mode: "copy"

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("${sample}.trimmed_merged.fastq.gz")

    script:
    """
    fastp \
        --in1 $read1 \
        --in2 $read2 \
        -p -c --merge \
        --merged_out ${sample}.trimmed_merged.fastq.gz \
        -h ${sample}.fastp_report.html \
        -j ${sample}.fastp_report.json \
        -w ${task.cpus} \
        -l 30
    """
}

process KRAKENUNIQ {
    tag "$sample:$db_name"
    publishDir "results/${sample}/krakenuniq", mode: "copy"

    input:
    tuple val(sample), path(trimmed) 
    val(db_name)
    val(db_dir)

    output:
    path "sequences_${sample}.krakenuniq_${db_name}"
    path "krakenuniq_${sample}.output_${db_name}"

    script:
    """
    krakenuniq \
        --db $db_dir \
        --fastq-input $trimmed \
        --threads ${task.cpus} \
        --output sequences_${sample}.krakenuniq_${db_name} \
        --report-file krakenuniq_${sample}.output_${db_name} \
        --preload-size ${task.memory.toGiga()}G
    """
}

workflow {
    // Make a channel from the samplesheet
    samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true, sep:'\t')
        .map { row ->
            tuple(row.sample, file(row.read1), file(row.read2))
        }

    // Run fastp
    trimmed_ch = samples_ch | FASTP

    // Define databases
    dbs = [
        [ "MAMMALIAN",   params.db_mammalian ],
        [ "OTHER",       params.db_other ],
        [ "PLANT",       params.db_plant ],
        [ "INVERTEBRATE",params.db_invertebrate ],
        [ "MICROBE",     params.db_microbe ]
    ]

    // Expand each sample across all databases
    trimmed_ch
        .combine(Channel.fromList(dbs))
        .map { sample, trimmed, db_info -> tuple(sample, trimmed, db_info[0], db_info[1]) }
        | KRAKENUNIQ
}

