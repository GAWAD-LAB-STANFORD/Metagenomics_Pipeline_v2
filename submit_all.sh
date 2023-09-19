#!/bin/bash
#
#SBATCH --job-name=submit_all
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=5:00:00
#SBATCH --partition=cgawad

PIPELINE_DIR="$( cd "$( dirname "$0" )" && pwd )"
PIPELINE_COMMAND="$@"
HELP="\
Purpose: \n\t\
    This pipeline is built to identify metagenomic species from pair-end fastq.gz files and remove human contamination \n\n\
Required arguments: -p/--project <arg> and either -f/--fastq_dir <arg> or -r/--results_dir <arg> \n\
Optional arguments: -s/--scratch_dir <arg>, --err_out_dir <arg>, --skip_scratch, -b/--run_dir <arg>, \n\t\
    --sample_sheet <arg>, --skip_identify, --only_identify, --identify <arg>, \n\t\
    --R1_suffix <arg>, --R2_suffix <arg>, --skip_trimming, --filter_rhesus, --kraken_db_types <arg>, --min_kraken_reads <arg>, \n\t\
    --subspecies, --blast_db_types <arg>, --num_alignments <arg>, --align_min <arg>, --slurm <arg> \n\
Defaults: \n\t\
    If no fastq_dir specified, uses results_dir \n\t\
    If no results_dir specified, makes new directory in fastq_dir \n\t\
    scratch_dir: /scratch/groups/cgawad/date_project_Scratch \n\t\
    sample_sheet: SampleSheet.csv \n\t\
    R1_suffix: _L001_R1_001.fastq.gz or _R1_001.fastq.gz or _R1.fastq.gz \n\t\
    R2_suffix: _L001_R2_001.fastq.gz or _R2_001.fastq.gz or _R2.fastq.gz \n\t\
    kraken_db_types: microbial \n\t\
    min_kraken_reads: 10 \n\t\
    blast_db_types: nt \n\t\
    num_alignments: 250 \n\t\
    align_min: 90 \n\n\
Run after demultiplexing and with fastq directory: \n\t\
    sh ${PIPELINE_DIR}/submit_all.sh --fastq_dir /oak/stanford/groups/cgawad/2020-01-01_Fastqs/ --project 2020-01-01_Project \n\n\
Run after demultiplexing and with results directory: \n\t\
    sh ${PIPELINE_DIR}/submit_all.sh --fastq_dir /oak/stanford/groups/cgawad/2020-01-01_Fastqs/ --results_dir /oak/stanford/groups/cgawad/2020-01-01_Results/ --project 2020-01-01_Project \n\n\
Run with demultiplexing and fastq directory: \n\t\
    sh ${PIPELINE_DIR}/submit_all.sh --run_dir /oak/stanford/groups/cgawad/Illumina_Data/MiniSeq/2020-01-01_BCLs --fastq_dir /oak/stanford/groups/cgawad/2020-01-01_Fastqs/ --project 2020-01-01_Project \n\n\
Run with demultiplexing and results directory: \n\t\
    sh ${PIPELINE_DIR}/submit_all.sh --run_dir /oak/stanford/groups/cgawad/Illumina_Data/MiniSeq/2020-01-01_BCLs --results_dir /oak/stanford/groups/cgawad/2020-01-01_Results/ --project 2020-01-01_Project \n\n\
Run with demultiplexing, fastq directory, and results directory: \n\t\
    sh ${PIPELINE_DIR}/submit_all.sh --run_dir /oak/stanford/groups/cgawad/Illumina_Data/MiniSeq/2020-01-01_BCLs --fastq_dir /oak/stanford/groups/cgawad/2020-01-01_Fastqs/ --results_dir /oak/stanford/groups/cgawad/2020-01-01_Results/ --project 2020-01-01_Project \n\n\
For more information, read the README.md"

# Reads in command line option arguments and assigns them to variables
SKIP_SCRATCH=0
SKIP_IDENTIFY=0
ONLY_IDENTIFY=0
SKIP_TRIMMOMATIC=0
FILTER_RHESUS=1
KRAKEN_DB_TYPES="microbial-plasmid-viral"
MIN_KRAKEN_READS="10"
SUBSPECIES=0
BLAST_DB_TYPES="nt-plasmid-viral"
NUM_ALIGNMENTS=250
ALIGN_MINIMUM=90
STEP=0
TEMP_ARRAY_START=0
DEPENDENCIES=()
while [ "$1" != "" ]; do
    case $1 in
        -h | --help )           echo -e $HELP
                                exit 0
                                ;;
        -f | --fastq_dir )      shift
                                FASTQ_DIR=$1
                                ;;
        -r | --results_dir )    shift
                                RESULTS_DIR=$1
                                ;;
        -p | --project )        shift
                                PROJECT=$1
                                ;;
        -d | --pipeline_dir )   shift
                                PIPELINE_DIR=$1
                                ;;
        -s | --scratch_dir )    shift
                                SCRATCH_DIR=$1
                                ;;
        --err_out_dir )         shift
                                STD_ERR_OUT_DIR=$1
                                ;;
        --skip_scratch )        SKIP_SCRATCH=1
                                ;;
        -b | --run_dir )        shift
                                RUN_DIR=$1
                                ;;
        --sample_sheet )        shift
                                SAMPLE_SHEET=$1
                                ;;
        --skip_identify )       SKIP_IDENTIFY=1
                                ;;
        --only_identify )       ONLY_IDENTIFY=1
                                ;;
        --identify )            shift
                                IDENTIFY=$1
                                ;;
        --R1_suffix )           shift
                                R1_SUFFIX=$1
                                ;;
        --R2_suffix )           shift
                                R2_SUFFIX=$1
                                ;;
        --skip_trimming )       SKIP_TRIMMOMATIC=1
                                ;;
        --filter_rhesus )       FILTER_RHESUS=1
                                ;;
        --kraken_db_types )     shift
                                KRAKEN_DB_TYPES=$1
                                ;;
        --min_kraken_reads )    shift
                                MIN_KRAKEN_READS=$1
                                ;;
        --subspecies )          SUBSPECIES=1
                                ;;
        --blast_db_types )      shift
                                BLAST_DB_TYPES=$1
                                ;;
        --num_alignments )      shift
                                NUM_ALIGNMENTS=$1
                                ;;
        --align_min )           shift
                                ALIGN_MINIMUM=$1
                                ;;
        --step1 )               STEP=1
                                ;;
        --step2 )               STEP=2
                                ;;
        --step3 )               STEP=3
                                ;;
        --temp_array_start )    shift
                                TEMP_ARRAY_START=$1
                                ;;
        --slurm )               shift
                                SLURM_OPTIONS=${@:1}
                                ;;
    esac
    shift
done

# Hardcoded paths and variables
TEMP_ARRAY_INCREMENT=1000
REFERENCE_DIR="/oak/stanford/groups/cgawad/Reference_Files/GATK_Resource_Bundle_hg38"
TOOLS_DIR="/oak/stanford/groups/cgawad/Sequencing_Analysis_Tools"
REF_FASTA="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
RHESUS_FASTA="/oak/stanford/groups/cgawad/Reference_Files/Macaca_mulattta_Rhesus_monkey_hg38/Macaca_mulatta_Rhesus_monkey_hg38.fasta"
SCRIPT_DIR="${PIPELINE_DIR}/scripts"
KRAKEN_DB_DIR_PREFIX="/oak/stanford/groups/cgawad/Reference_Files/Kraken2_Fatfree_Databases/kraken2-fatfree-"
NCBI_DB_DIR_PREFIX="/oak/stanford/groups/cgawad/Reference_Files/NCBI_RefSeq_Databases/ncbi_database_"
NCBI_ANNOTATIONS_DIR="/oak/stanford/groups/cgawad/Reference_Files/NCBI_Annotations"

# Ensure we have the required variables set and set other variables
if ([ -z $FASTQ_DIR ] && [ -z $RESULTS_DIR ]) || [ -z $PROJECT ] || [ -z $PIPELINE_DIR ]; then
    echo "Variables not supplied correctly. Use -h/--help options for assistance. Ending program..."
    exit 1
fi
if [ -z $FASTQ_DIR ]; then
    FASTQ_DIR="$RESULTS_DIR"
elif [ -z $RESULTS_DIR ]; then
    RESULTS_DIR="${FASTQ_DIR}/$(date '+%Y-%m-%d')_${PROJECT}_Results"
fi
if [ ! -z $SCRATCH_DIR ] && [ $SKIP_SCRATCH -eq 1 ]; then
    echo "Variables not supplied correctly. Cannot skip scratch while being provided scratch_dir for use. Exiting with code 1"
    exit 1
fi
if [ -z $SCRATCH_DIR ] && [ $SKIP_SCRATCH -eq 0 ]; then
    SCRATCH_DIR="/scratch/groups/cgawad/$(date '+%Y-%m-%d')_${PROJECT}_Scratch"
fi
if [ -z $SCRATCH_DIR ] && [ $SKIP_SCRATCH -eq 1 ]; then
    SCRATCH_DIR="$RESULTS_DIR"
fi
if [ -z $STD_ERR_OUT_DIR ]; then
    STD_ERR_OUT_DIR="${RESULTS_DIR}/std_err_out_files"
fi
# Make directories if they don't exist
if [ ! -d $FASTQ_DIR ]; then
    mkdir $FASTQ_DIR
fi
if [ ! -d $RESULTS_DIR ]; then
    mkdir $RESULTS_DIR
fi
if [ ! -d $SCRATCH_DIR ]; then
    mkdir $SCRATCH_DIR
fi
if [ ! -d $STD_ERR_OUT_DIR ]; then
    mkdir $STD_ERR_OUT_DIR
fi
OPTIONS=( "-f $FASTQ_DIR -r $RESULTS_DIR -d $PIPELINE_DIR -p $PROJECT -s $SCRATCH_DIR --err_out_dir $STD_ERR_OUT_DIR " )
if [ ! -z $FASTQ_DIR ]; then
    OPTIONS+=( "-f $FASTQ_DIR" )
fi
if [ ! -z $RUN_DIR ] && [ $ONLY_IDENTIFY -eq 1 ]; then
    echo "Variables not supplied correctly. Cannot perform demultiplexing while only identifying data. Exiting with code 1"
    exit 1
fi
if [ ! -z $RUN_DIR ]; then
    SAMPLE_SHEET="${RUN_DIR}/SampleSheet.csv"
elif [ ! -z $SAMPLE_SHEET ]; then
    echo "Variables not supplied correctly. Please specify a run diretory for demultiplexing with --run_dir. Exiting with code 1"
    exit 1
fi
if [ ! -z $RUN_DIR ] && [ ! -z $SAMPLE_SHEET ]; then
    if [ ! -f $SAMPLE_SHEET ]; then
        echo "Sample sheet $SAMPLE_SHEET not found. Exiting with code 1"
        exit 1
    fi
    OPTIONS+=( "--run_dir $RUN_DIR --sample_sheet $SAMPLE_SHEET" )
fi
if [ $SKIP_IDENTIFY -eq 1 ] && [ $ONLY_IDENTIFY -eq 1 ]; then
    echo "Variables not supplied correctly. Please specify either --skip_identify or --only_identify, not both. Exiting with code 1"
    exit 1
elif [ $SKIP_IDENTIFY -eq 1 ]; then
    OPTIONS+=( "--skip_identify" )
elif [ $ONLY_IDENTIFY -eq 1 ]; then
    OPTIONS+=( "--only_identify" )
    if [ $STEP -eq 0 ]; then
        STEP=2
    fi
fi
if [ ! -z $IDENTIFY ]; then
    OPTIONS+=( "--identify $IDENTIFY" )
else
    IDENTIFY=$PROJECT
fi
if [ ! -z $R1_SUFFIX ]; then
    OPTIONS+=( "--R1_suffix $R1_SUFFIX" )
fi
if [ ! -z $R2_SUFFIX ]; then
    OPTIONS+=( "--R2_suffix $R2_SUFFIX" )
fi
if [ $SKIP_TRIMMOMATIC -eq 1 ]; then
    OPTIONS+=( "--skip_trimming" )
fi

if [ $FILTER_RHESUS -eq 1 ]; then
    OPTIONS+=( "--filter_rhesus" )
    REF_FASTA_STRING="${REF_FASTA}:${RHESUS_FASTA}"
    REF_NAME_STRING="human:rhesus"
else
    REF_FASTA_STRING="${REF_FASTA}"
    REF_NAME_STRING="human"
fi
if [ "$KRAKEN_DB_TYPES" != "microbial" ]; then
    OPTIONS+=( "--kraken_db_types $KRAKEN_DB_TYPES" )
fi
if [ $STEP -eq 0 ] && [ -z $RUN_DIR ]; then
    STEP=1
fi
if [ "$BLAST_DB_TYPES" != "nt" ]; then
    OPTIONS+=( "--blast_db_types $BLAST_DB_TYPES" )
fi
if [ $NUM_ALIGNMENTS -ne 250 ]; then
    OPTIONS+=( "--num_alignments $NUM_ALIGNMENTS" )
fi
if [ $ALIGN_MINIMUM -ne 90 ]; then
    OPTIONS+=( "--align_min $ALIGN_MINIMUM" )
fi


TEMP_PIPELINE_DIR="$( cd "$( dirname "$0" )" && pwd )"
if [ $ONLY_IDENTIFY -eq 1 ]; then
    PIPELINE_STATUS=${STD_ERR_OUT_DIR}/${IDENTIFY}_pipeline_status.txt
else
    PIPELINE_STATUS=${STD_ERR_OUT_DIR}/${PROJECT}_pipeline_status.txt
fi
cd $SCRATCH_DIR
if [ "$TEMP_PIPELINE_DIR" = "$PIPELINE_DIR" ]; then
    echo -e "\nSTART: $(date)\nMetagenomics Pipeline v2\n\n$PIPELINE_COMMAND\n\nProject: $PROJECT\nResults dir: $RESULTS_DIR\nFastq dir: $FASTQ_DIR\nScratch dir: $SCRATCH_DIR\nErr out dir: $STD_ERR_OUT_DIR" >> $PIPELINE_STATUS
    # Optional variable definitions
    if [ $SKIP_SCRATCH -eq 0 ]; then
        echo "Default: Scratch dir is different from Results dir" >> $PIPELINE_STATUS
    else
        echo "Option: Scratch dir is the same as Results dir" >> $PIPELINE_STATUS
    fi
    if [ $SKIP_IDENTIFY -eq 1 ]; then
        echo "Option: Skip identification of data - will only process the fastqs and skip BAM processing" >> $PIPELINE_STATUS
    fi
    if [ $ONLY_IDENTIFY -eq 1 ]; then
        echo "Option: Only identification of data - will only process already constructed BAMs" >> $PIPELINE_STATUS
    fi
    if [ "$IDENTIFY" != "$PROJECT" ]; then
        echo "Option: Identify different from Project: $IDENTIFY" >> $PIPELINE_STATUS
    fi
    if [ $SKIP_TRIMMOMATIC -eq 1 ]; then
        echo "Option: Skip trimming - will not run trimmomatic" >> $PIPELINE_STATUS
    fi
    if [ $FILTER_RHESUS -eq 1 ]; then
        echo "Option: Filter rhesus - will remove reads that align to macaca mulatta rhesus monkey" >> $PIPELINE_STATUS
    fi
    if [ "$KRAKEN_DB_TYPES" = "microbial" ]; then
        echo "Default: Kraken db types: microbial" >> $PIPELINE_STATUS
    else
        echo "Option: Kraken db types: $KRAKEN_DB_TYPES" >> $PIPELINE_STATUS
    fi
    if [ "$MIN_KRAKEN_READS" = "10" ]; then
        echo "Default: Minimum kraken reads: 10" >> $PIPELINE_STATUS
    else
        echo "Option: Minimum kraken reads: $MIN_KRAKEN_READS" >> $PIPELINE_STATUS
    fi
    if [ $SUBSPECIES -eq 1 ]; then
        echo "Option: Subspecies - will include subspecies (kraken S1) in analysis" >> $PIPELINE_STATUS
    fi
    if [ "$BLAST_DB_TYPES" = "nt" ]; then
        echo "Default: Blast db types: nt" >> $PIPELINE_STATUS
    else
        echo "Option: Blast db types: $BLAST_DB_TYPES" >> $PIPELINE_STATUS
    fi
    if [ $NUM_ALIGNMENTS -eq 250 ]; then
        echo "Default: Num alignments: 250" >> $PIPELINE_STATUS
    else
        echo "Option: Num alignments: $NUM_ALIGNMENTS" >> $PIPELINE_STATUS
    fi
    if [ $ALIGN_MINIMUM -eq 90 ]; then
        echo "Default: Align minimum: 90" >> $PIPELINE_STATUS
    else
        echo "Option: Align minimum: $ALIGN_MINIMUM" >> $PIPELINE_STATUS
    fi
    echo " " >> $PIPELINE_STATUS
fi


if [ ! -z $SLURM_OPTIONS ]; then
    echo "Slurm option used - entire pipeline run will be queued with user parameters" >> $PIPELINE_STATUS
    sbatch -J $PROJECT ${SLURM_OPTIONS[@]} \
        -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
        ${PIPELINE_DIR}/submit_all.sh ${OPTIONS[@]}
    exit 0
fi


if [ $STEP -ne 0 ] && [ $ONLY_IDENTIFY -eq 0 ]; then
    if [ -z $R1_SUFFIX ] || [ -z $R2_SUFFIX ]; then
        R1_SUFFIX="_L001_R1_001.fastq.gz"
        R2_SUFFIX="_L001_R2_001.fastq.gz"
        if [ $(find ${FASTQ_DIR} -maxdepth 1 -name "*${R1_SUFFIX}" | wc -l) -eq 0 ]; then
            R1_SUFFIX="_R1_001.fastq.gz"
            R2_SUFFIX="_R2_001.fastq.gz"
        fi
        if [ $(find ${FASTQ_DIR} -maxdepth 1 -name "*${R1_SUFFIX}" | wc -l) -eq 0 ]; then
            R1_SUFFIX="_R1.fastq.gz"
            R2_SUFFIX="_R2.fastq.gz"
        fi
    fi
    SAMPLE_ARRAY=( $(find ${FASTQ_DIR} -maxdepth 1 -name "*${R1_SUFFIX}" -exec basename {} \; | \
        grep -v "Undetermined" | sed "s/${R1_SUFFIX}//") )
    if [ ${#SAMPLE_ARRAY[@]} -eq 0 ]; then
        echo "No fastq.gz files found in the fastq directory. Exiting with code 1" >> $PIPELINE_STATUS
        echo "END: $(date)" >> $PIPELINE_STATUS
        exit 1
    fi
    if [ $STEP -eq 1 ]; then
        echo -e "Number of samples: ${#SAMPLE_ARRAY[@]}\nSamples: ${SAMPLE_ARRAY[@]}" >> $PIPELINE_STATUS
    fi
fi


if [ $STEP -eq 0 ]; then
    echo "### Demultiplexing ### - START: $(date)" >> $PIPELINE_STATUS
    echo -e "Run dir: $RUN_DIR\nSample sheet: $SAMPLE_SHEET" >> $PIPELINE_STATUS
    echo -e "\nsbatch --parsable -e ${STD_ERR_OUT_DIR}/%A_%x.err -o ${STD_ERR_OUT_DIR}/%A_%x.out \
        ${SCRIPT_DIR}/0_demultiplexer.sh --run_dir $RUN_DIR --sample_sheet $SAMPLE_SHEET --fastq_dir $FASTQ_DIR \
        --pipeline_status $PIPELINE_STATUS\n" >> $PIPELINE_STATUS
    DEPENDENCIES+=( $(sbatch --parsable -e ${STD_ERR_OUT_DIR}/%A_%x.err -o ${STD_ERR_OUT_DIR}/%A_%x.out \
        ${SCRIPT_DIR}/0_demultiplexer.sh --run_dir $RUN_DIR --sample_sheet $SAMPLE_SHEET --fastq_dir $FASTQ_DIR \
        --pipeline_status $PIPELINE_STATUS) )
    echo -e "\nsbatch --dependency=afterok:${DEPENDENCIES[0]} -J $PROJECT \
        -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
        ${PIPELINE_DIR}/submit_all.sh --step1 ${OPTIONS[@]}\n" >> $PIPELINE_STATUS
    sbatch --dependency=afterok:${DEPENDENCIES[0]} -J $PROJECT \
        -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
        ${PIPELINE_DIR}/submit_all.sh --step1 ${OPTIONS[@]}
elif [ $STEP -eq 1 ]; then
    if [ $TEMP_ARRAY_START -eq 0 ]; then
        echo "### Processing fastq samples ### - START: $(date)" >> $PIPELINE_STATUS
        JOB_COUNT=${#SAMPLE_ARRAY[@]}
        echo "Process sample jobs to run: $JOB_COUNT" >> $PIPELINE_STATUS
        TEMP_ARRAY_START=1
    fi

    TEMP_SAMPLE_ARRAY=( ${SAMPLE_ARRAY[@]:$(($TEMP_ARRAY_START - 1)):$TEMP_ARRAY_INCREMENT} )
    TEMP_JOB_COUNT=${#TEMP_SAMPLE_ARRAY[@]}
    echo -e "$(date)\nSubmitting $TEMP_JOB_COUNT jobs for samples $TEMP_ARRAY_START to $(($TEMP_ARRAY_START + ${#TEMP_SAMPLE_ARRAY[@]} - 1))" >> $PIPELINE_STATUS
    TEMP_SAMPLES_STRING=$( IFS=$':'; echo "${TEMP_SAMPLE_ARRAY[*]}" )
    echo -e "\nsbatch --parsable -e $STD_ERR_OUT_DIR/%A_%a_%x.err -o $STD_ERR_OUT_DIR/%A_%a_%x.out \
        --array=1-${TEMP_JOB_COUNT} ${SCRIPT_DIR}/1_process_fastqs.sh \
        $FASTQ_DIR $SCRATCH_DIR $R1_SUFFIX $R2_SUFFIX $SKIP_TRIMMOMATIC \
        $REF_FASTA_STRING $REF_NAME_STRING $TOOLS_DIR $TEMP_SAMPLES_STRING\n" >> $PIPELINE_STATUS
    DEPENDENCIES+=( $(sbatch --parsable -e $STD_ERR_OUT_DIR/%A_%a_%x.err -o $STD_ERR_OUT_DIR/%A_%a_%x.out \
        --array=1-${TEMP_JOB_COUNT} ${SCRIPT_DIR}/1_process_fastqs.sh \
        $FASTQ_DIR $SCRATCH_DIR $R1_SUFFIX $R2_SUFFIX $SKIP_TRIMMOMATIC \
        $REF_FASTA_STRING $REF_NAME_STRING $TOOLS_DIR $TEMP_SAMPLES_STRING) )
    TEMP_ARRAY_START=$(($TEMP_ARRAY_START + $TEMP_ARRAY_INCREMENT))
    echo -e "$(date)\nIncrement: $TEMP_ARRAY_INCREMENT\nNew start: $TEMP_ARRAY_START" >> $PIPELINE_STATUS
    
    if [ $TEMP_ARRAY_START -le ${#SAMPLE_ARRAY[@]} ]; then
        echo -e "\nsbatch --dependency=afterok:${DEPENDENCIES[0]} -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step1 --temp_array_start $TEMP_ARRAY_START ${OPTIONS[@]}\n"  >> $PIPELINE_STATUS
        sbatch --dependency=afterok:${DEPENDENCIES[0]} -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step1 --temp_array_start $TEMP_ARRAY_START ${OPTIONS[@]}
    else
        echo -e "\nsbatch --dependency=afterok:${DEPENDENCIES[0]} -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step2 ${OPTIONS[@]}\n"  >> $PIPELINE_STATUS
        sbatch --dependency=afterok:${DEPENDENCIES[0]} -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step2 ${OPTIONS[@]}
    fi
elif [ $STEP -eq 2 ] && [ $ONLY_IDENTIFY -eq 0 ]; then
    if [ $TEMP_ARRAY_START -eq 0 ]; then
        SAMPLE_COUNT=1
        for SAMPLE in ${SAMPLE_ARRAY[@]}; do
            if [ ! -f ${SAMPLE}_ref_filtered.bam ]; then
                echo "Sample number $SAMPLE_COUNT - ${SAMPLE}_ref_filtered.bam not found. Exiting with code 1" >> $PIPELINE_STATUS
                echo "END: $(date)" >> $PIPELINE_STATUS
                exit 1
            fi
            SAMPLE_COUNT=$((SAMPLE_COUNT+1))
        done
        rm ${STD_ERR_OUT_DIR}/*1_process_fastqs.out ${STD_ERR_OUT_DIR}/*1_process_fastqs.err
        echo "### Processing fastq samples ### - END: $(date)" >> $PIPELINE_STATUS
    fi

    if [ $SKIP_IDENTIFY -eq 1 ]; then
        echo "Ending without identfication of data" >> $PIPELINE_STATUS
        echo "END: $(date)" >> $PIPELINE_STATUS
        exit 0
    fi
fi


if [ $STEP -eq 2 ] || [ $STEP -eq 3 ]; then
    SAMPLE_ARRAY=( $(ls *_ref_filtered.bam | sed "s/_ref_filtered.bam//") )
    if [ ${#SAMPLE_ARRAY[@]} -eq 0 ]; then
        echo "No filtered BAM files found in the results directory. Exiting with code 1" >> $PIPELINE_STATUS
        echo "END: $(date)" >> $PIPELINE_STATUS
        exit 1
    fi
fi


if [ $STEP -eq 2 ]; then
    if [ $TEMP_ARRAY_START -eq 0 ]; then
        echo "### Processing ref filtered samples ### - START: $(date)" >> $PIPELINE_STATUS
        JOB_COUNT=${#SAMPLE_ARRAY[@]}
        echo "Process sample jobs to run: $JOB_COUNT" >> $PIPELINE_STATUS
        TEMP_ARRAY_START=1
    fi

    TEMP_SAMPLE_ARRAY=( ${SAMPLE_ARRAY[@]:$(($TEMP_ARRAY_START - 1)):$TEMP_ARRAY_INCREMENT} )
    TEMP_JOB_COUNT=${#TEMP_SAMPLE_ARRAY[@]}
    echo -e "$(date)\nSubmitting $TEMP_JOB_COUNT jobs for samples $TEMP_ARRAY_START to $(($TEMP_ARRAY_START + ${#TEMP_SAMPLE_ARRAY[@]} - 1))" >> $PIPELINE_STATUS
    TEMP_SAMPLES_STRING=$( IFS=$':'; echo "${TEMP_SAMPLE_ARRAY[*]}" )
    echo -e "\nsbatch --parsable -e $STD_ERR_OUT_DIR/%A_%a_%x.err -o $STD_ERR_OUT_DIR/%A_%a_%x.out \
        --array=1-${TEMP_JOB_COUNT} ${SCRIPT_DIR}/2_process_sample.sh \
        $SCRATCH_DIR $R1_SUFFIX $R2_SUFFIX $KRAKEN_DB_TYPES $KRAKEN_DB_DIR_PREFIX $BLAST_DB_TYPES \
        $NCBI_DB_DIR_PREFIX $NUM_ALIGNMENTS $ALIGN_MINIMUM $TOOLS_DIR $SCRIPT_DIR $PROJECT \
        $MIN_KRAKEN_READS $SUBSPECIES $TEMP_SAMPLES_STRING\n" >> $PIPELINE_STATUS
    DEPENDENCIES+=( $(sbatch --parsable -e $STD_ERR_OUT_DIR/%A_%a_%x.err -o $STD_ERR_OUT_DIR/%A_%a_%x.out \
        --array=1-${TEMP_JOB_COUNT} ${SCRIPT_DIR}/2_process_sample.sh \
        $SCRATCH_DIR $R1_SUFFIX $R2_SUFFIX $KRAKEN_DB_TYPES $KRAKEN_DB_DIR_PREFIX $BLAST_DB_TYPES \
        $NCBI_DB_DIR_PREFIX $NUM_ALIGNMENTS $ALIGN_MINIMUM $TOOLS_DIR $SCRIPT_DIR $PROJECT \
        $MIN_KRAKEN_READS $SUBSPECIES $TEMP_SAMPLES_STRING) )
    TEMP_ARRAY_START=$(($TEMP_ARRAY_START + $TEMP_ARRAY_INCREMENT))
    echo -e "New start: $TEMP_ARRAY_START\nIncrement: $TEMP_ARRAY_INCREMENT" >> $PIPELINE_STATUS
    
    if [ $TEMP_ARRAY_START -le ${#FASTQ_ARRAY[@]} ]; then
        echo -e "\nsbatch --dependency=afterany:$( IFS=$':'; echo "${DEPENDENCIES[*]}" ) -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step2 --temp_array_start $TEMP_ARRAY_START ${OPTIONS[@]}\n" >> $PIPELINE_STATUS
        sbatch --dependency=afterany:$( IFS=$':'; echo "${DEPENDENCIES[*]}" ) -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step2 --temp_array_start $TEMP_ARRAY_START ${OPTIONS[@]}
    else
        echo -e "\nsbatch --dependency=afterany:$( IFS=$':'; echo "${DEPENDENCIES[*]}" ) -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step3 ${OPTIONS[@]}\n" >> $PIPELINE_STATUS
        sbatch --dependency=afterany:$( IFS=$':'; echo "${DEPENDENCIES[*]}" ) -J $PROJECT \
            -e ${STD_ERR_OUT_DIR}/%A_submit_all_%x.err -o ${STD_ERR_OUT_DIR}/%A_submit_all_%x.out \
            ${PIPELINE_DIR}/submit_all.sh --step3 ${OPTIONS[@]}
    fi
elif [ $STEP -eq 3 ]; then
    KRAKEN_DB_TYPE_ARRAY=( $(echo $KRAKEN_DB_TYPES | sed 's/-/ /g') )
    SAMPLE_COUNT=1
    for SAMPLE in ${SAMPLE_ARRAY[@]}; do
        if [ ! -f ${PROJECT}.${SAMPLE}_${KRAKEN_DB_TYPE_ARRAY[0]}_kraken_report.tsv ]; then
            echo -e "\tSample number $SAMPLE_COUNT - ${PROJECT}.${SAMPLE}.${KRAKEN_DB_TYPE_ARRAY[0]}_kraken_report.tsv file not found" >> $PIPELINE_STATUS
        fi
        SAMPLE_COUNT=$((SAMPLE_COUNT+1))
    done
    RESULTS_COUNT=$(ls ${PROJECT}.*_kraken_report.tsv | wc -l)
    MAX_RESULTS=$(echo ${#SAMPLE_ARRAY[@]} ${#KRAKEN_DB_TYPE_ARRAY[@]} | awk '{ print $1 * $2 }')
    if [ $RESULTS_COUNT -eq 0 ]; then
        echo "No results found. Exiting with code 1" >> $PIPELINE_STATUS
        echo "END: $(date)" >> $PIPELINE_STATUS
        exit 1
    else
        echo "$RESULTS_COUNT results out of a possible $MAX_RESULTS maximum" >> $PIPELINE_STATUS
        rm ${STD_ERR_OUT_DIR}/*2_process_sample.out ${STD_ERR_OUT_DIR}/*2_process_sample.err
    fi
    echo "### Processing ref filtered samples ### - END: $(date)" >> $PIPELINE_STATUS
    
    
    echo "### Consolidating files ### - START: $(date)" >> $PIPELINE_STATUS
    if [ ! -z $FASTQ_DIR ]; then
        REF_NAME_ARRAY=( $(echo $REF_NAME_STRING | sed 's/:/ /g') )
        for REF_NAME in ${REF_NAME_ARRAY[@]}; do
            ALIGNMENT_METRICS_FILENAMES=( $(ls *_${REF_NAME}_alignment_metrics.tsv) )
            echo -e sample"\t"$(head -n 7 ${ALIGNMENT_METRICS_FILENAMES[0]} | tail -n 1) | sed 's/ /\t/g' > ${PROJECT}.${REF_NAME}_alignment_metrics_final.tsv
            for i in ${ALIGNMENT_METRICS_FILENAMES[@]}; do 
                SAMPLE=$(echo $i | sed "s/_${REF_NAME}_alignment_metrics.tsv//")
                R1=$(head -n 8 $i | tail -n 1)
                R2=$(head -n 9 $i | tail -n 1)
                PAIR=$(head -n 10 $i | tail -n 1)
                echo -e "$SAMPLE\t$R1\n$SAMPLE\t$R2\n$SAMPLE\t$PAIR"
            done | sed 's/ /\t/g' >> ${PROJECT}.${REF_NAME}_alignment_metrics_final.tsv
            echo "Merged $REF_NAME alignment metrics" >> $PIPELINE_STATUS
            
            rm ${ALIGNMENT_METRICS_FILENAMES[@]}
        done
    fi
    
    KRAKEN_DB_TYPE_ARRAY=( $(echo $KRAKEN_DB_TYPES | sed 's/-/ /g') )
    for DB_TYPE in ${KRAKEN_DB_TYPE_ARRAY[@]}; do
        KRAKEN_REPORT_FILENAMES=( $(ls ${PROJECT}.*_${DB_TYPE}_kraken_report.tsv) )
        echo -e "sample\tpercent_fragments_covered\tfragments_covered\tfragments_assigned\trank_code\ttaxid\tsciname" > \
            ${PROJECT}.${DB_TYPE}_kraken_reports_final.tsv
        for i in ${KRAKEN_REPORT_FILENAMES[@]}; do
            SAMPLE=$(echo $i | sed "s/_${DB_TYPE}_kraken_report.tsv//" | sed "s/${PROJECT}.//")
            cat $i | sed 's/^ \+/'${SAMPLE}'\t/' >> ${PROJECT}.${DB_TYPE}_kraken_reports_final.tsv
        done
        
        KRAKEN_BLAST_FILENAMES=( $(ls ${PROJECT}.*_${DB_TYPE}_kraken_blast.tsv) )
        head -n 1 ${KRAKEN_BLAST_FILENAMES[0]} > ${PROJECT}.${DB_TYPE}_kraken_blast_final.tsv
        for i in ${KRAKEN_BLAST_FILENAMES[@]}; do
            tail -n +2 $i >> ${PROJECT}.${DB_TYPE}_kraken_blast_final.tsv
        done
        
        rm ${KRAKEN_REPORT_FILENAMES[@]} ${KRAKEN_BLAST_FILENAMES[@]}
    done
    echo "### Consolidating files ### - END: $(date)" >> $PIPELINE_STATUS
    
    
    # echo "### Converting Kraken reports to TSV ### - START: $(date)" >> $PIPELINE_STATUS
    # KRAKEN_DB_TYPE_ARRAY=( $(echo $KRAKEN_DB_TYPES | sed 's/-/ /g') )
    # SAMPLES_STRING=$( IFS=$':'; echo "${SAMPLE_ARRAY[*]}" )
    # echo "Samples string: $SAMPLES_STRING"
    # for DB_TYPE in ${KRAKEN_DB_TYPE_ARRAY[@]}; do
    #     python3 ${SCRIPT_DIR}/kraken_report_to_jtree.py ${PROJECT}.${DB_TYPE}.kraken_reports.tsv \
    #         $IDENTIFY .${DB_TYPE}.kraken_jtree.json $SAMPLES_STRING
    # done
    # echo "### Converting Kraken reports to TSV ### - END: $(date)" >> $PIPELINE_STATUS
    
    
    # echo "### Processing consolidated results and making figures ### - START: $(date)" >> $PIPELINE_STATUS
    # ml R/4.2.0
    # export R_LIBS="/home/groups/cgawad/R_LIBS"
    # Rscript ${SCRIPT_DIR}/analyze_and_plot_results.R \
    #     --project $PROJECT --identify $IDENTIFY \
    #     --sample_read_count_filename ${PROJECT}.sample_read_counts.tsv \
    #     --summed_read_targets_filename ${PROJECT}.summed_read_targets.tsv \
    #     --contig_read_targets_filename ${PROJECT}.contig_read_targets.tsv \
    #     --contig_data_filename ${PROJECT}.contig_data.tsv \
    #     --kraken_db_types $KRAKEN_DB_TYPES --kraken_jtree_suffix ".kraken_jtree.json" \
    #     --blast_db_types $BLAST_DB_TYPES --blast_results_suffix ".blast_results.tsv" \
    #     --ncbi_annotations_dir $NCBI_ANNOTATIONS_DIR --contig_alignment_percent_min $CONTIG_ALIGN_MINIMUM
    # echo "### Processing consolidated results and making figures ### - END: $(date)" >> $PIPELINE_STATUS
    
    if [ "$SCRATCH_DIR" != "$RESULTS_DIR" ]; then
        echo "### Moving results from scratch dir to results dir ### - START: $(date)"
        rsync -ar $SCRATCH_DIR $RESULTS_DIR
        echo "### Moving results from scratch dir to results dir ### - END: $(date)"
    fi
    echo "END: $(date)" >> $PIPELINE_STATUS
fi
