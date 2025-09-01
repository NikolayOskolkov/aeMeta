# aeMeta: ancient environmental metagenomic workflow

This is a computational workflow for fast profiling ancient environmental metagenomic data. The workflow accepts paired-end data in FASTQ-format and outputs a table of abundances of eukaryotic organisms. The workflow uses pre-built KrakenUniq databases and Bowtie2 indexes constructed for NCBI RefSeq plant, vertebrate, mammals, and invertabrate chromosome-level genome assemblies.

The workflow was developed by Nikolay Oskolkov, Lund University, Sweden.

If you use the workflow for your research, please cite our manuscript:

    Nikolay Oskolkov,
    Profiling ancient environmental metagenomes using aeMeta 
    https://www.biorxiv.org/content/10.1101/2025.03.19.644176v1, 
    https://doi.org/10.1101/2025.03.19.644176

Please note that in this gitub reporsitory, we provide a small subset of microbial pseudo-reads for demonstration purposes, the full dataset is available at the SciLifeLab Figshare https://doi.org/10.17044/scilifelab.28491956.

Questions regarding the dataset should be sent to nikolay.oskolkov@scilifelab.se

## Quick start

Please clone this repository and install the workflow tools as follows:

    git clone https://github.com/NikolayOskolkov/aeMeta
    cd aeMeta
    conda env create -f environment.yaml
    conda activate aeMeta

Then you can specify the workflow input files and parameters in the `nextflow.config` and run it using Nextflow:

    nextflow run main.nf

The Nextflow implementation is preferred for scalability and reproducibility purposes. Please place your reference genomes (fasta-files) to be screened for exogenous regions in the `data` folder. An example of the config-file, `nextflow.config`, can look like this:

    params {
        input_dir = "data"                                             // folder with multiple reference genomes (fasta-files)
        type_of_pseudo_reads = "GTDB"                                  // type of pseudo-reads to be used for screening the input reference genome, can be "GTDB", "RefSeq" or "human"
        threads = 4                                                    // number of available threads
        input_pseudo_reads = "GTDB_sliced_seqs_sliding_window.fna.gz"  // name of pre-computed file with pseudo-reads, can be "GTDB_sliced_seqs_sliding_window.fna.gz", "RefSeq_sliced_seqs_sliding_window.fna.gz" or "human_sliced_seqs_sliding_window.fna.gz"
        n_allowed_multimappers = 10                                    // number of multi-mapping pseudo-reads allowed by Bowtie2, do not change this default number unless you know what you are doing
    }

Please modify it to adjust for the number of available threads in your computational environment and the type of analysis, i.e. detecting microbial-like or human-like sequeneces in the reference genome, you would like to perform.
