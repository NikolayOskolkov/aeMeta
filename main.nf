#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.samplesheet = ""

params.db_mammalian    = "DBDIR_VERTEBRATE_MAMMALIAN"
params.db_other        = "DBDIR_VERTEBRATE_OTHER"
params.db_plant        = "DBDIR_PLANT"
params.db_invertebrate = "DBDIR_INVERTEBRATE"
params.db_microbe      = "DBDIR_MICROBE"

// -------------------------
// FASTP process
// -------------------------
process FASTP {
    tag "$sample"
    publishDir "results/${sample}/fastp", mode: "copy"

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), 
          path("${sample}.trimmed_merged.fastq.gz"), 
          path("${sample}.fastp_report.html"), 
          path("${sample}.fastp_report.json")

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

// -------------------------
// KRAKENUNIQ process
// -------------------------
process KRAKENUNIQ {
    tag "$sample:$db_name"
    publishDir "results/${sample}/krakenuniq", mode: "copy"

    input:
    tuple val(sample), path(trimmed), val(db_name), val(db_dir)

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

// -------------------------
// KRAKENUNIQ_ABUNDANCE process
// -------------------------
process KRAKENUNIQ_ABUNDANCE {
    tag "abundance_calc"
    publishDir "results/abundance", mode: "copy"

    input:
    path kraken_dirs  // all directories collected

    output:
    path "abundance_matrix_eukaryotes.txt"

    script:
    """
    Rscript /home/nikolay-oskolkov/WABI/HPGG/aeMeta/scripts/krakenuniq_abundance.R
    """
}

// -------------------------
// Workflow
// -------------------------
workflow {

    // 1. Read samplesheet
    samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true, sep:'\t')
        .map { row -> tuple(row.sample, file(row.read1), file(row.read2)) }

    // 2. Run FASTP
    trimmed_ch = samples_ch | FASTP

    // 3. Define databases
    db_list = [
        ["MAMMALIAN",    params.db_mammalian],
        ["OTHER",        params.db_other],
        ["PLANT",        params.db_plant],
        ["INVERTEBRATE", params.db_invertebrate],
        ["MICROBE",      params.db_microbe]
    ]

    // 4. Prepare KRAKENUNIQ inputs
    kraken_input_ch = trimmed_ch
        .map { sample, trimmed, html, json -> tuple(sample, trimmed) }
        .flatMap { sample, trimmed ->
            db_list.collect { db_name, db_dir ->
                tuple(sample, trimmed, db_name, db_dir)
            }
        }

    // 5. Run KRAKENUNIQ
    kraken_input_ch | KRAKENUNIQ

    // 6. Collect all KRAKENUNIQ output directories AFTER all runs complete
    kraken_dirs_ch = Channel
        .fromPath("results/*/krakenuniq")
        .collect()

    // 7. Run abundance calculation
    kraken_dirs_ch | KRAKENUNIQ_ABUNDANCE
}

