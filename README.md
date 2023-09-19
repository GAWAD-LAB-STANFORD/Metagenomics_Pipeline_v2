# Guide to Metagenomics_Pipeline

- [Purpose](#purpose)
- [How To Run](#how-to-run)
- [What It Does Exactly](#what-it-does-exactly)
- [Resources](#resources)
- [TODO and Notes](#todo-and-notes)

## Purpose
- This pipeline is built to identify metagenomic species from pair-end fastq.gz files and remove human contamination

## How To Run
- All scripts are controlled by the master script, pipeline_control.sh
- To make job submission easier, use the submit_all.sh script to submit the master script
    - For help, use the *-h* or *--help* option like:

    ```bash
    sh submit_all.sh --help
    ```

- You can use ~ if your file or folder is in your home directory, and you can exclude a path if your file or folder is in the results directory, but otherwise use absolute paths
- You can have '/' or nothing at the end of a directory path, either is fine:
    - -r /home/groups/cgawad/results/
    - -r /home/groups/cgawad/results

### submit_all.sh
- **Run this script in order to run the entire pipeline**
- Required arguments: -p/--project >arg< and either -f/--fastq_dir >arg< or -r/--results_dir >arg<
    - Specify where your fastq.gz files are located (or will be put after you opt for the auto-demultiplexing) using *-f* or *--fastq_dir* and/or specify where to output the results using *-r* or *--results_dir*
        - If you do not specify a fastq directory, the program will assume the fastq.gz files are in the results directory you specified, and will end the program if no fastq.gz files are found
        - If you do not specify a results directory, the program will make a new folder with the current date in the name within the fastq directory 
    - Specify the project name for the final resulting VCF that will be made using *-p* or *--project*
- Optional arguments: -s/--scratch_dir >arg<, --err_out_dir >arg<, --skip_scratch, -b/--run_dir >arg<, --sample_sheet >arg<, --skip_identify, --only_identify, --identify >arg<, --R1_suffix >arg<, --R2_suffix >arg<, --skip_trimming, --rna, --filter_rhesus, --kraken_db_types >arg<, --min_kraken_reads >arg<, --subspecies, --blast_db_types >arg<, --num_alignments >arg<, --align_min >arg<, --slurm >arg<
    - You can specify a directory to perform all intermediate steps in with *-s* or *--scratch_dir*
    - You can specify a directory to output the standard error and out print statements of all jobs to using *--err_out_dir*
    - If you want to skip having the pipeline run intermediate steps in scratch, use *--skip_scratch*
    - You can have the script demultiplex your BCL files into fastq.gz files by specifying a run folder using *-b* or *--run_dir*. The program will look for a sample sheet called SampleSheet.csv in the first level within the run_dir or you can specify a different sample sheet with *--sample_sheet*. The program will make the fastq directory if it does not exist and tell you the sizes of undeteremined vs fully demultiplexed reads
    - If you only want to process the fastqs, build the contigs, and run Kraken2, add the *--skip_identify* option and only the first phase of the pipeline will be run
    - If you want to BLAST your data, filter the results, and make the final figures using already existing contig files, you can use the *--only_identify*
    - If you want to name the results from the 2nd half the pipeline differently from the 1st half, which would be helpful in the case where you want to perform numerous different and simultaneous identification and filtering runs on the same data, add the *--identify* option with a different argument from *--project*
    - If your read 1 and read 2 fastq.gz files differentiate themselves by some pattern other than _L001_R1_001.fastq.gz and _L001_R2_001.fastq.gz or _R1_001.fastq.gz and _R2_001.fastq.gz or _R1.fastq.gz and R2.fastq.gz, use *--R1_suffix* and *--R2_suffix* options to let the pipeline know
    - If you don't want trimmomatic to run, add the *--skip_trimming* option
    - If you want to process RNA instead of DNA data, add the *--rna* option and the program will align the sequences using STAR instead of BWA
    - If you want to also filter reads from Rhesus monkey, add the *--filter_rhesus* option
    - If you want to use different kraken2 databases from the default, specify which ones you want to use with *--kraken_db_types*
    - If you want to filter kraken reads by some other number, use the *--min_kraken_reads* option
    - If you want to include subspecies along with species in the analysis, add the *--subspecies* option
    - If you want to use different BLAST databases from the default, specify which ones you want to use with *--blast_db_types*
    - If you want to only have returned a maximum number of BLAST results, choose your maximum number using *--num_alignments*
    - You can change the default minimum alignment for BLAST results from 0.9 (90%) to a chosen number using *--align_min*
    - Besides the already implemented job name and standard error and output print statements, you can specify additional slurm commands for the pipeline job following the use of the *--slurm* option. If you use this option, **make sure it is the last one you use**
        - A useful example would be setting a future time to run the job and asking for email notifications like so:

        ```
        ... --slurm --begin=now+12hours --mail-type=ALL 
        ```

- Get some example submissions by using *-h* or *--help* options like:

    ```bash
    sh submit_all.sh --help
    ```

## What It Does Exactly

![Pipeline Graphic](pipeline_graphic.png)

### submit_all.sh
- **0_demultiplexer.sh** - Demultiplexing will be done if specified
- **1_process_fastqs.sh** - For each sample, a job will be run to do the following:
    - Trimmomatic trimming
    - BWN ALN align all reads to human genome
    - Separate non-human reads
    - Collect trimmomatic, read, and alignment metrics
- **2_process_sample.sh**
    - Kraken2 run on non-human reads using specified kraken2 libraries to identify species
    - SPAdes de novo assembly of contigs using non-human reads and correlated Kraken2 species
    - BLAST+ (blastn) of species contigs using specified blast libraries to classify contigs
    - Parse blast results
    - Consolidate species results across kraken2 and blast libraries
- Consolidate alignment metrics and kraken2 results
- Create json jtree files from kraken2 report outputs
- Process consolidated metrics and create figures


## Resources
- How do I create a database for BLASTing?
    - Download NCBI database fasta files using a command like this:

    ```bash
    rsync --copy-links --recursive --times --verbose --progress rsync://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria*.genomic.fna.gz ./
    ```

    - Combine all those fasta files into an uncompressed single fasta file
    - Make a BLAST database from that combined single fasta file using a command like this:

    ```bash
    makeblastdb -in combined_bacteria_sequences.fasta -dbtype nucl -out bacteria
    ```

    - If you use NCBI's own BLAST database instead of creating your own (their databases are more condensed and take up less space), then download and update using a command like this:

    ```bash
    perl update_blastdb.pl --decompress nt
    ```

- How do I get annotations for my sequences?
    - Download NCBI database gbff annotation files using a command like this:

    ```bash
    rsync --copy-links --recursive --times --verbose --progress rsync://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria*.genomic.gbff.gz ./
    ```

    - Combine all those gbff files into an uncompressed single gbff file (or to do this faster with parallelization, do this step for each gbff file individually and combine the parsed annotations at the very end), then run a script like the one I made (parse_genbank_annotations.py) with your modifications

    ```bash
    python3 parse_genbank_annotations.py uncompressed_bacteria.gbff parsed_bacteria_annotations.tsv
    ```

## TODO and Notes