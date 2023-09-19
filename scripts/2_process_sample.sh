#!/bin/bash
#
#SBATCH --job-name=2_process_sample
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=4-00:00:00
#SBATCH --partition=cgawad

START_TIME=$(date +%s)
SCRATCH_DIR=$1
FASTQ_DIR=$2
R1_SUFFIX=$3
R2_SUFFIX=$4
KRAKEN_DB_TYPE_ARRAY=( $(echo $5 | sed 's/-/ /g') )
KRAKEN_DB_DIR_PREFIX=$6
BLAST_DB_TYPE_ARRAY=( $(echo $7 | sed 's/-/ /g') )
NCBI_DB_DIR_PREFIX=$8
NUM_ALIGNMENTS=$9
ALIGN_MINIMUM=${10}
TOOLS_DIR=${11}
SCRIPT_DIR=${12}
PROJECT=${13}
MIN_KRAKEN_READS=${14}
SUBSPECIES=${15}
SAMPLE_ARRAY=( $(echo ${16} | sed 's/:/ /g') )
SAMPLE=${SAMPLE_ARRAY[$(( $SLURM_ARRAY_TASK_ID - 1 ))]}

echo -e "START: $(date)\nMetagenomics pipeline\nResults dir: $SCRATCH_DIR\nSample: $SAMPLE"
cd $SCRATCH_DIR

ml java perl R/4.2.0 python/3.6.1 py-pandas/0.23.0_py36 py-numpy/1.14.3_py36
ml biology bwa samtools gatk
export R_LIBS="/home/groups/cgawad/R_LIBS"
export PATH=${TOOLS_DIR}/kraken2-2.0.8-beta:$PATH

samtools view ${SAMPLE}_ref_filtered.bam > ${PROJECT}.${SAMPLE}_ref_filtered.sam
# samtools view ${SAMPLE}_no_rhesus.bam > ${SAMPLE}_ref_filtered.sam
for DB_TYPE in ${KRAKEN_DB_TYPE_ARRAY[@]}; do
    echo "### Identifying matches between unaligned reads and kraken2 $DB_TYPE database ### - START: $(date)"
    kraken2 --db ${KRAKEN_DB_DIR_PREFIX}${DB_TYPE} --threads 4 --fastq-input --paired --gzip-compressed \
        --output ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_vs_ref_filtered.tsv  \
        --report ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_report.tsv \
        ${FASTQ_DIR}/${SAMPLE}${R1_SUFFIX} ${FASTQ_DIR}/${SAMPLE}${R2_SUFFIX}
    echo "### Identifying matches between unaligned reads and kraken2 $DB_TYPE database ### - END: $(date)"
    
    echo "### Filtering out species from kraken2 $DB_TYPE results ### - START: $(date)"
    ### Bacteria and Fungi analysis before January 11, 2021 used direct kraken reads >= 50 without subspecies (only S, not S1)
    if [ $SUBSPECIES -eq 1 ]; then
        cat ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_report.tsv | \
            awk -v pat=$MIN_KRAKEN_READS '{ if ($3 >= pat && ($4 == "S" || $4 == "S1")) { print } }' | \
            tr -s '  ' ' ' | sed 's/^[ ]*//' > ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_species.tsv
    else
        cat ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_report.tsv | \
            awk -v pat=$MIN_KRAKEN_READS '{ if ($3 >= pat && $4 == "S") { print } }' | \
            tr -s '  ' ' ' | sed 's/^[ ]*//' > ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_species.tsv
    fi
    echo "### Filtering out species from kraken2 $DB_TYPE results ### - END: $(date)"
    
    echo "### Building species contigs from $DB_TYPE results ### - START: $(date)"
    SPECIES_COUNT=$(cat ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_species.tsv | wc -l)
    for ((SPECIES_NUM = 1 ; SPECIES_NUM <= $SPECIES_COUNT ; SPECIES_NUM++)); do
        SPECIES_ID=$(sed "${SPECIES_NUM}q;d" ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_species.tsv | cut -f 5)
        SPECIES=$(sed "${SPECIES_NUM}q;d" ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_species.tsv | cut -f 6 | sed "s/^[ ]*//")
        echo -e "START: $(date)\nSpecies ID $SPECIES_ID - Species number $SPECIES_NUM of $SPECIES_COUNT"
        awk -v pat=$SPECIES_ID '$3 == pat' ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_vs_ref_filtered.tsv | \
            cut -f 2 > ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_read_headers.tsv
        echo -e "\tSpecies read headers filtered"
        
        grep -wFf ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_read_headers.tsv \
            ${PROJECT}.${SAMPLE}_ref_filtered.sam > ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}.sam
        echo -e "\tSpecies reads filtered"
        
        gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx64g" SamToFastq \
            -I ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}.sam --VALIDATION_STRINGENCY SILENT \
            -F ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}${R1_SUFFIX} \
            -F2 ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}${R2_SUFFIX}
        echo -e "\tSpecies converted from reads to fastqs"
        
        ${TOOLS_DIR}/BBTools/bbmap_38.87/dedupe.sh \
            in=${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}${R1_SUFFIX} \
            out=${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R1.fasta s=3
        awk 'NR % 3 == 1' ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R1.fasta | \
            sed 's/>//g' | sed "s/.$/2/g" > ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R1_header
        zcat ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}${R2_SUFFIX} | \
            awk "NR % 4 == 1 || NR % 4 == 2" | sed 's/@/>/g' \
            > ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R2_raw.fasta
        grep -A 1 -wFf  ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R1_header \
            ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R2_raw.fasta \
            >${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R2.fasta
        cat ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R1.fasta \
            ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R2.fasta | sed "s/--//g" | awk 'NF' \
            > ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}.fasta
        echo -e "\tSpecies converted from fastqs to fastas"
        
        COUNT_BLAST_DB_TYPE=1
        NUM_BLAST_DB_TYPES=${#BLAST_DB_TYPE_ARRAY[@]}
        for BLAST_DB_TYPE in ${BLAST_DB_TYPE_ARRAY[@]}; do
            echo -e "\tSTART: $(date)\nBlast DB type $BLAST_DB_TYPE - Blast DB number $COUNT_BLAST_DB_TYPE of $NUM_BLAST_DB_TYPES"
            export BLASTDB=${NCBI_DB_DIR_PREFIX}${BLAST_DB_TYPE}
            ${TOOLS_DIR}/ncbi-blast-2.10.0+/bin/blastn -db $BLAST_DB_TYPE -num_alignments ${NUM_ALIGNMENTS} \
                -num_threads 2 -outfmt 15 -query ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}.fasta \
                -out ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_${BLAST_DB_TYPE}_blast.json
            echo -e "\t\tSpecies blasted for blast db type $BLAST_DB_TYPE"
            
            python3 ${SCRIPT_DIR}/blast_json_to_tsv.py \
                -i ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_${BLAST_DB_TYPE}_blast.json \
                -o ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_${BLAST_DB_TYPE}_blast.tsv \
                -b $BLAST_DB_TYPE -s $SAMPLE -k $DB_TYPE -d $SPECIES_ID
            echo -e "\t\tSpecies $BLAST_DB_TYPE blast results converted from json to tsv"
            
            Rscript ${SCRIPT_DIR}/process_blast_tsv.R \
                ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_${BLAST_DB_TYPE}_blast.tsv $ALIGN_MINIMUM
            if [ $(cat ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_${BLAST_DB_TYPE}_blast.tsv | wc -l) -le 1 ]; then
                rm ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_${BLAST_DB_TYPE}_blast.tsv
            fi
            echo -e "\t\tSpecies $BLAST_DB_TYPE blast results processed"
             
            rm ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_${BLAST_DB_TYPE}_blast.json
            COUNT_BLAST_DB_TYPE=$((COUNT_BLAST_DB_TYPE+1))
        done
        
        rm ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_read_headers.tsv
        rm ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}.sam
        rm ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}${R1_SUFFIX} 
        rm ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}${R2_SUFFIX}
        rm ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}.fasta
        rm ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R2_raw.fasta
        rm ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R1.fasta
        rm ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R2.fasta
        rm ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R1_header
    done
    echo "### Building species contigs from $DB_TYPE results ### - END: $(date)"
    rm ${PROJECT}.${SAMPLE}_ref_filtered.sam
    
    echo "### Consolidating final species files for $DB_TYPE results ### - START: $(date)"
    BLAST_FILENAMES=( $(ls ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_*_blast.tsv) )
    head -n 1 ${BLAST_FILENAMES[0]} > ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_blast.tsv
    for i in ${BLAST_FILENAMES[@]}; do
        tail -n +2 $i >> ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_blast.tsv
    done
    if [ $(cat ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_blast.tsv | wc -l) -le 1 ]; then
        rm ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_blast.tsv
    else
        python3 ${SCRIPT_DIR}/merge_kraken_blast.py \
            -b ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_blast.tsv \
            -k ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_report.tsv
    fi
    echo "${#BLAST_FILENAMES[@]} species blast json files consolidated"
    
    rm ${BLAST_FILENAMES[@]}
    rm ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_species.tsv
    rm ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_vs_ref_filtered.tsv
    echo "### Consolidating final species files for $DB_TYPE results ### - END: $(date)"
done

if [ ! -f ${PROJECT}.${SAMPLE}_${KRAKEN_DB_TYPE_ARRAY[0]}_kraken_report.tsv ]; then
    echo "Final file ${PROJECT}.${SAMPLE}_${KRAKEN_DB_TYPE_ARRAY[0]}_kraken_report.tsv not found. Exiting with code 1"
    exit 1
fi
echo -e "END: $(date)\nRuntime: $(($(date +%s)-$START_TIME)) seconds"
