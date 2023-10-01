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
TOOLS_DIR=$5
PROJECT=$6
KRAKEN_DB_TYPE_ARRAY=( $(echo $7 | sed 's/-/ /g') )
KRAKEN_DB_DIR_PREFIX=$8
MIN_KRAKEN_READS=$9
SUBSPECIES=${10}
BLAST_CONTIGS=${11}
BLAST_DB_TYPE_ARRAY=( $(echo ${12} | sed 's/-/ /g') )
NCBI_DB_DIR_PREFIX=${13}
NUM_ALIGNMENTS=${14}
SCRIPT_DIR=${15}
ALIGN_MINIMUM=${16}
SAMPLE_ARRAY=( $(echo ${17} | sed 's/:/ /g') )
SAMPLE=${SAMPLE_ARRAY[$(( $SLURM_ARRAY_TASK_ID - 1 ))]}

echo -e "START: $(date)\nMetagenomics pipeline v2\nResults dir: $SCRATCH_DIR\nSample: $SAMPLE"
cd $SCRATCH_DIR

ml java perl R/4.2.0 python/3.6.1 py-pandas/0.23.0_py36 py-numpy/1.14.3_py36
ml biology bwa samtools gatk
export R_LIBS="/home/groups/cgawad/R_LIBS"
export PATH=${TOOLS_DIR}/kraken2-2.0.8-beta:$PATH
export PATH=${TOOLS_DIR}/mmseqs/bin/:$PATH

samtools view ${SAMPLE}_ref_filtered.bam > temp_${SAMPLE}_ref_filtered.sam
if [ $BLAST_CONTIGS -eq 1 ]; then
    READS=$(samtools view ${SAMPLE}_ref_filtered.bam | cut -f 3 | grep "chr" | wc -l)
    echo -e "sample\ttotal_reads\tspecies_id\tkraken_db\tcontig_aligned\tcontig_unaligned" > ${PROJECT}.${SAMPLE}_read_targets_counts.tsv
    echo -e "read_count\tsample\tspecies_id\tkraken_db\ttarget" > ${PROJECT}.${SAMPLE}_kraken_contig_read_targets.tsv
fi
COUNT_KRAKEN_DB_TYPE=1
NUM_KRAKEN_DB_TYPES=${#KRAKEN_DB_TYPE_ARRAY[@]}
for DB_TYPE in ${KRAKEN_DB_TYPE_ARRAY[@]}; do
    echo "$COUNT_KRAKEN_DB_TYPE of $NUM_KRAKEN_DB_TYPES Kraken DB types - Type: $DB_TYPE - START: $(date)"
    echo -e "\t### Identifying matches between unaligned reads and kraken2 $DB_TYPE database ### - START: $(date)"
    kraken2 --db ${KRAKEN_DB_DIR_PREFIX}${DB_TYPE} --threads 4 --fastq-input --paired --gzip-compressed \
        --output temp_${SAMPLE}_${DB_TYPE}_kraken_vs_ref_filtered.tsv  \
        --report ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_report.tsv \
        ${FASTQ_DIR}/${SAMPLE}${R1_SUFFIX} ${FASTQ_DIR}/${SAMPLE}${R2_SUFFIX}
    echo -e "\t### Identifying matches between unaligned reads and kraken2 $DB_TYPE database ### - END: $(date)"
    
    
    echo -e "\t### Filtering out species from kraken2 $DB_TYPE results ### - START: $(date)"
    ### Bacteria and Fungi analysis before January 11, 2021 used direct kraken reads >= 50 without subspecies (only S, not S1)
    if [ $SUBSPECIES -eq 1 ]; then
        cat ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_report.tsv | \
            awk -v pat=$MIN_KRAKEN_READS '{ if ($3 >= pat && ($4 == "S" || $4 == "S1")) { print } }' | \
            tr -s '  ' ' ' | sed 's/^[ ]*//' > temp_${SAMPLE}_${DB_TYPE}_kraken_species.tsv
    else
        cat ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_report.tsv | \
            awk -v pat=$MIN_KRAKEN_READS '{ if ($3 >= pat && $4 == "S") { print } }' | \
            tr -s '  ' ' ' | sed 's/^[ ]*//' > temp_${SAMPLE}_${DB_TYPE}_kraken_species.tsv
    fi
    echo -e "\t### Filtering out species from kraken2 $DB_TYPE results ### - END: $(date)"
    
    
    echo -e "\t### Comparing kraken species from $DB_TYPE with BLAST ### - START: $(date)"
    SPECIES_COUNT=$(cat temp_${SAMPLE}_${DB_TYPE}_kraken_species.tsv | wc -l)
    for ((SPECIES_NUM = 1 ; SPECIES_NUM <= $SPECIES_COUNT ; SPECIES_NUM++)); do
        SPECIES_ID=$(sed "${SPECIES_NUM}q;d" temp_${SAMPLE}_${DB_TYPE}_kraken_species.tsv | cut -f 5)
        SPECIES=$(sed "${SPECIES_NUM}q;d" temp_${SAMPLE}_${DB_TYPE}_kraken_species.tsv | cut -f 6 | sed "s/^[ ]*//")
        echo -e "\t$SPECIES_NUM of $SPECIES_COUNT species - Species ID: $SPECIES_ID - START: $(date)"
        awk -v pat=$SPECIES_ID '$3 == pat' temp_${SAMPLE}_${DB_TYPE}_kraken_vs_ref_filtered.tsv | \
            cut -f 2 > temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_read_headers.tsv
        echo -e "\t\tSpecies read headers filtered"
        
        grep -wFf temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_read_headers.tsv \
            temp_${SAMPLE}_ref_filtered.sam > temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}.sam
        echo -e "\t\tSpecies reads filtered"
        
        gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx64g" SamToFastq \
            -I temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}.sam --VALIDATION_STRINGENCY SILENT \
            -F temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}${R1_SUFFIX} \
            -F2 temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}${R2_SUFFIX}
        echo -e "\t\tSpecies converted from reads to fastqs"
        
        if [ $BLAST_CONTIGS -eq 1 ]; then
            mkdir spades_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}
            python3 /oak/stanford/groups/cgawad/Sequencing_Analysis_Tools/SPAdes-3.14.0-Linux/bin/spades.py \
                -t 4 -m 64 -1 temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}${R1_SUFFIX} -2 temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}${R2_SUFFIX} \
                -o spades_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}
            if [ ! -f spades_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}/contigs.fasta ]; then
                echo -e "\t\tWARNING: No contigs made, skipping rest of this kraken species"
                continue
            fi
            mv spades_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}/contigs.fasta temp_${SAMPLE}.${SPECIES_ID}_${DB_TYPE}_kraken_contigs.fasta
            echo -e "\t\tSpecies contigs built"
            
            bwa index temp_${SAMPLE}.${SPECIES_ID}_${DB_TYPE}_kraken_contigs.fasta
            bwa mem -M temp_${SAMPLE}.${SPECIES_ID}_${DB_TYPE}_kraken_contigs.fasta \
                ${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}${R1_SUFFIX} ${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}${R2_SUFFIX} | \
                samtools view -b - | samtools sort -o ${SAMPLE}.${SPECIES_ID}_${DB_TYPE}_kraken_contig_aligned.bam -
            echo -e "\t\tAligned reads to contigs"
            
            gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx64g" CollectAlignmentSummaryMetrics \
                -R temp_${SAMPLE}.${SPECIES_ID}_${DB_TYPE}_kraken_contigs.fasta -I temp_${SAMPLE}_${SPECIES_ID}_${DB_TYPE}_kraken_contig_aligned.bam -O temp_${SAMPLE}_${SPECIES_ID}_${DB_TYPE}_contig_alignment_metrics.tsv
            echo -e "\t\tCollected read-to-contig alignment metrics"
            
            samtools view temp_${SAMPLE}_${SPECIES_ID}_${DB_TYPE}_kraken_contig_aligned.bam | cut -f 3 | \
                grep "NODE" > temp_${SAMPLE}_${SPECIES_ID}_${DB_TYPE}_kraken_contig_read_targets.txt
            printf "$SAMPLE\t$SPECIES_ID\t$DB_TYPE\n%0.s" $(seq $(cat temp_${SAMPLE}_${SPECIES_ID}_${DB_TYPE}_kraken_contig_read_targets.txt | wc -l)) | \
                paste - temp_${SAMPLE}-${SPECIES_ID}_${DB_TYPE}_kraken_contig_read_targets.txt | uniq -c | \
                sed 's/^[[:space:]]*//' | tr -s ' ' '\t'  >> temp_${SAMPLE}_kraken_contig_read_targets.tsv
            echo -e "\t\tCollected read-to-contig target data"
            
            CONTIG_ALIGNED_READS=$(samtools view temp_${SAMPLE}_${SPECIES_ID}_${DB_TYPE}_kraken_contig_aligned.bam | cut -f 3 | grep "NODE" | wc -l)
            CONTIG_UNALIGNED_READS=$(samtools view temp_${SAMPLE}_${SPECIES_ID}_${DB_TYPE}_kraken_contig_aligned.bam | cut -f 3 | grep -v "NODE" | wc -l)
            echo -e "$SAMPLE\t$READS\t$SPECIES_ID\t$DB_TYPE\t$CONTIG_ALIGNED_READS\t$CONTIG_UNALIGNED_READS" >> ${PROJECT}.${SAMPLE}_read_targets_counts.tsv
            echo -e "\t\tCollected read-to-contig target counts"
        else
            ${TOOLS_DIR}/BBTools/bbmap_38.87/dedupe.sh \
                in=temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}${R1_SUFFIX} \
                out=temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R1.fasta s=3
            awk 'NR % 3 == 1' temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R1.fasta | \
                sed 's/>//g' | sed "s/.$/2/g" > temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R1_header
            zcat temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}${R2_SUFFIX} | \
                awk "NR % 4 == 1 || NR % 4 == 2" | sed 's/@/>/g' \
                > temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R2_raw.fasta
            grep -A 1 -wFf  temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R1_header \
                temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R2_raw.fasta \
                > temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R2.fasta
            cat temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R1.fasta \
                temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R2.fasta | sed "s/--//g" | awk 'NF' \
                > temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}.fasta
            echo -e "\t\tSpecies converted from fastqs to fastas"
        fi
        
        COUNT_BLAST_DB_TYPE=1
        NUM_BLAST_DB_TYPES=${#BLAST_DB_TYPE_ARRAY[@]}
        for BLAST_DB_TYPE in ${BLAST_DB_TYPE_ARRAY[@]}; do
            echo -e "\t\t$COUNT_BLAST_DB_TYPE of $NUM_BLAST_DB_TYPES BLAST DB types - Blast DB type: $BLAST_DB_TYPE - START: $(date)"
            export BLASTDB=${NCBI_DB_DIR_PREFIX}${BLAST_DB_TYPE}
            
            if [ $BLAST_CONTIGS -eq 0 ]; then
                ${TOOLS_DIR}/ncbi-blast-2.10.0+/bin/blastn -db $BLAST_DB_TYPE -num_alignments ${NUM_ALIGNMENTS} \
                    -num_threads 2 -outfmt 15 -query temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}.fasta \
                    -out ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_${BLAST_DB_TYPE}_blast.json
            else
                ${TOOLS_DIR}/ncbi-blast-2.10.0+/bin/blastn -db $BLAST_DB_TYPE -num_alignments ${NUM_ALIGNMENTS} \
                    -num_threads 2 -outfmt 15 -query temp_${SAMPLE}.${SPECIES_ID}_${DB_TYPE}_kraken_contigs.fasta \
                    -out ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_${BLAST_DB_TYPE}_blast.json
            fi
            echo -e "\t\t\tSpecies blasted"
            
            python3 ${SCRIPT_DIR}/blast_json_to_tsv.py \
                -i ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_${BLAST_DB_TYPE}_blast.json \
                -o ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_${BLAST_DB_TYPE}_blast.tsv \
                -b $BLAST_DB_TYPE -s $SAMPLE -k $DB_TYPE -d $SPECIES_ID
            echo -e "\t\t\tSpecies blast results converted from json to tsv"
            
            Rscript ${SCRIPT_DIR}/process_blast_tsv.R \
                ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_${BLAST_DB_TYPE}_blast.tsv $ALIGN_MINIMUM
            if [ $(cat ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_${BLAST_DB_TYPE}_blast.tsv | wc -l) -le 1 ]; then
                rm ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_${BLAST_DB_TYPE}_blast.tsv
            fi
            echo -e "\t\t\tSpecies blast results processed"
             
            rm ${PROJECT}.${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_${BLAST_DB_TYPE}_blast.json
            COUNT_BLAST_DB_TYPE=$((COUNT_BLAST_DB_TYPE+1))
        done
        
        rm temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_read_headers.tsv
        rm temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}.sam
        rm temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}${R1_SUFFIX} 
        rm temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}${R2_SUFFIX}
        if [ $BLAST_CONTIGS -eq 0 ]; then
            rm temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}.fasta
            rm temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R2_raw.fasta
            rm temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R1.fasta
            rm temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R2.fasta
            rm temp_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}_R1_header
        else
            rm -r spades_${SAMPLE}_${DB_TYPE}_kraken_${SPECIES_ID}
            rm temp_${SAMPLE}_${SPECIES_ID}_${DB_TYPE}_kraken_contig_aligned.bam
            rm temp_${SAMPLE}_${SPECIES_ID}_${DB_TYPE}_kraken_contig_aligned.bam.bai
            rm temp_${SAMPLE}_${SPECIES_ID}_${DB_TYPE}_kraken_temp_contig_read_targets.txt
        fi
    done
    rm temp_${SAMPLE}_${DB_TYPE}_kraken_species.tsv
    rm temp_${SAMPLE}_${DB_TYPE}_kraken_vs_ref_filtered.tsv
    COUNT_KRAKEN_DB_TYPE=$((COUNT_KRAKEN_DB_TYPE+1))
    echo -e "\t### Comparing kraken species from $DB_TYPE with BLAST ### - END: $(date)"
    
    
    echo -e "\t### Consolidating final species files for $DB_TYPE results ### - START: $(date)"
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
    echo -e "\t### Consolidating final species files for $DB_TYPE results ### - END: $(date)"
done
rm temp_${SAMPLE}_ref_filtered.sam


if [ $BLAST_CONTIGS -eq 1 ]; then
    echo "### Consolidating contig alignment metrics ### - START: $(date)"
    CONTIG_ALIGNMENT_METRICS_FILENAMES=( $(ls temp_${SAMPLE}_*_contig_alignment_metrics.tsv) )
    echo -e "sample\tspecies_id\tkraken_db\t$(sed -n '7p' ${CONTIG_ALIGNMENT_METRICS_FILENAMES[0]})" > ${PROJECT}.${SAMPLE}_contig_alignment_metrics.tsv
    for i in ${CONTIG_ALIGNMENT_METRICS_FILENAMES[@]}; do
        SPECIES_ID=$(echo $i | cut -d '_' -f 3)
        DB_TYPE=$(echo $i | cut -d '_' -f 4)
        R1=$(head -n '8p' $i)
        R2=$(head -n '9p' $i)
        PAIR=$(head -n '10p' $i)
        echo -e "$SAMPLE\t$SPECIES_ID\t$DB_TYPE\t$R1\n$SAMPLE\t$SPECIES_ID\t$DB_TYPE\t$R2\n$SAMPLE\t$SPECIES_ID\t$DB_TYPE\t$PAIR\n" >> ${PROJECT}.${SAMPLE}_contig_alignment_metrics.tsv
    done
    rm ${CONTIG_ALIGNMENT_METRICS_FILENAMES[@]}
    echo "### Consolidating contig alignment metrics ### - END: $(date)"


    echo "### Consolidating contig data ### - START: $(date)"
    echo -e "sample\tcontig" > ${PROJECT}.contig_data.tsv
    CONTIGS_FILENAMES=( $(ls temp_${SAMPLE}_*_kraken_contigs.fasta) )
    for i in ${CONTIGS_FILENAMES[@]}; do
        SPECIES_ID=$(echo $i | cut -d '_' -f 3)
        DB_TYPE=$(echo $i | cut -d '_' -f 4
        grep ">" $i | xargs -i echo -e "${SAMPLE}\t${SPECIES_ID}\t${DB_TYPE}\t{}" >> ${PROJECT}.${SAMPLE}_contig_data.tsv
    done
    echo "Merged contigs" >> $PIPELINE_STATUS
    # rm ${CONTIGS_FILENAMES[@]}
    echo "### Consolidating contig data ### - END: $(date)"
fi


if [ ! -f ${PROJECT}.${SAMPLE}_${KRAKEN_DB_TYPE_ARRAY[0]}_kraken_report.tsv ]; then
    echo "Final file ${PROJECT}.${SAMPLE}_${KRAKEN_DB_TYPE_ARRAY[0]}_kraken_report.tsv not found. Exiting with code 1"
    exit 1
fi
echo -e "END: $(date)\nRuntime: $(($(date +%s)-$START_TIME)) seconds"
