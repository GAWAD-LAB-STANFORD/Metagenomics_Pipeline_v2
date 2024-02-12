#!/bin/bash
#
#SBATCH --job-name=demultiplexer
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=6:00:00
#SBATCH --partition=cgawad

START_TIME=$(date +%s)
SCRIPT_COMMAND="$@"
while [ "$1" != "" ]; do
    case $1 in
        --run_dir )             shift
                                RUN_DIR=$1
                                ;;
        --sample_sheet )        shift
                                SAMPLE_SHEET=$1
                                ;;
        --fastq_dir )           shift
                                FASTQ_DIR=$1
                                ;;
        --pipeline_status )     shift
                                PIPELINE_STATUS=$1
                                ;;
    esac
    shift
done

if [ -z $RUN_DIR ] || [ -z $SAMPLE_SHEET ] || [ -z $FASTQ_DIR ] || [ -z $PIPELINE_STATUS ]; then
    echo "Variables not supplied correctly. Use -h/--help options for assistance. Ending program..."
    exit 1
fi

echo -e "START: $(date)\nMetagenomics Pipeline v2\nScript command: $SCRIPT_COMMAND"

ml biology bcl2fastq
bcl2fastq --runfolder-dir $RUN_DIR --sample-sheet $SAMPLE_SHEET --output-dir $FASTQ_DIR

echo -e "END: $(date)\nRuntime: $(($(date +%s)-$START_TIME)) seconds"
if [ $(find $FASTQ_DIR -name "*.fastq.gz" | wc -l) -eq 0 ]; then
    echo "No fastq.gz files found in fastq dir. Exiting with code 1" >> $PIPELINE_STATUS
    exit 1
fi
find $FASTQ_DIR -name "*.fastq.gz" -exec mv {} ${FASTQ_DIR}/ \;
BCL_SIZE=$(du -sh $RUN_DIR | cut -f 1)
UND_SIZE=$(du -shc ${FASTQ_DIR}/Undetermined*.fastq.gz | tail -n 1 | cut -f 1)
FASTQ_SIZE=$(ls ${FASTQ_DIR}/*.fastq.gz | grep -v "Undetermined\|extracted" | xargs du -shc | tail -n 1 | cut -f 1)
echo -e "$BCL_SIZE run dir produced $UND_SIZE of undetermined fastq.gz and $FASTQ_SIZE of determined fastq.gz" >> $PIPELINE_STATUS
echo "### Demultiplexing ### - END: $(date)" >> $PIPELINE_STATUS