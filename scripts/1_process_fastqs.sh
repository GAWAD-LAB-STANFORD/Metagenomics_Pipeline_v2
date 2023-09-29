#!/bin/bash
#
#SBATCH --job-name=1_process_fastqs
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=4-00:00:00
#SBATCH --partition=cgawad

START_TIME=$(date +%s)
FASTQ_DIR=$1
SCRATCH_DIR=$2
R1_SUFFIX=$3
R2_SUFFIX=$4
SKIP_TRIMMOMATIC=$5
TOOLS_DIR=$6
RNA=$7
REF_FASTA_ARRAY=( $(echo $8 | sed 's/:/ /g') )
REF_NAME_ARRAY=( $(echo $9 | sed 's/:/ /g') )
SAMPLE_ARRAY=( $(echo ${10} | sed 's/:/ /g') )
SAMPLE=${SAMPLE_ARRAY[$(( $SLURM_ARRAY_TASK_ID - 1 ))]}

echo -e "START: $(date)\nMetagenomics pipeline v2\nFastq dir: $FASTQ_DIR\nResults dir: $SCRATCH_DIR\nSample: $SAMPLE"
cd $SCRATCH_DIR

ml java/11.0.11 perl R/4.2.0 python/3.6.1 py-pandas/0.23.0_py36 py-numpy/1.14.3_py36
ml biology bwa samtools gatk
export R_LIBS="/home/groups/cgawad/R_LIBS"

R1_FASTQ=${FASTQ_DIR}/${SAMPLE}${R1_SUFFIX}
R2_FASTQ=${FASTQ_DIR}/${SAMPLE}${R2_SUFFIX}

if [ $SKIP_TRIMMOMATIC -eq 0 ]; then
    UNTRIMMED_R1_FASTQ=$R1_FASTQ
    UNTRIMMED_R2_FASTQ=$R2_FASTQ
    R1_FASTQ=$(echo ${SAMPLE}${R1_SUFFIX} | sed "s/_R1/_R1_trimmed/")
    R2_FASTQ=$(echo ${SAMPLE}${R2_SUFFIX} | sed "s/_R2/_R2_trimmed/")
    UNPAIRED_R1_FASTQ=$(echo ${SAMPLE}${R1_SUFFIX} | sed "s/_R1/_R1_trimmed_unpaired/")
    UNPAIRED_R2_FASTQ=$(echo ${SAMPLE}${R2_SUFFIX} | sed "s/_R2/_R2_trimmed_unpaired/")
    
    echo "### Trimming fastqs ### - START: $(date)"
    java -jar ${TOOLS_DIR}/Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 -trimlog \
        ${SAMPLE}_trimmomatic_log.txt ${UNTRIMMED_R1_FASTQ} ${UNTRIMMED_R2_FASTQ} \
        ${R1_FASTQ} ${UNPAIRED_R1_FASTQ} \
        ${R2_FASTQ} ${UNPAIRED_R2_FASTQ} \
        ILLUMINACLIP:${TOOLS_DIR}/Trimmomatic-0.35/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
        LEADING:3 TRAILING:3 MINLEN:36
    echo "### Trimming fastqs ### - END: $(date)"
fi

echo "### Counting fastq read counts ### - START: $(date)"
READ_COUNT=$(echo $(zcat $R1_FASTQ | wc -l ) \
    $(zcat $R2_FASTQ | wc -l) | awk '{ print ($1 + $2) / 4 }' )
echo -e "sample\tread_count" > ${SAMPLE}.read_counts.tsv
echo -e "$SAMPLE\t$READ_COUNT" >> ${SAMPLE}.read_counts.tsv
echo "### Counting fastq read counts ### - END: $(date)"

PREV_R1_FASTQ=$R1_FASTQ
PREV_R2_FASTQ=$R2_FASTQ
for ((REF_INDEX = 0 ; REF_INDEX < ${#REF_FASTA_ARRAY[@]} ; REF_INDEX++)); do
    REF_FASTA=${REF_FASTA_ARRAY[$REF_INDEX]}
    REF_NAME=${REF_NAME_ARRAY[$REF_INDEX]}

    if [ $RNA -eq 1 ]; then
        echo "### Aligning RNA fastqs to $REF_NAME ### - START: $(date)"
        UNZIPPED_R1_FASTQ=$(basename $PREV_R1_FASTQ | sed "s/.gz//")
        UNZIPPED_R2_FASTQ=$(basename $PREV_R1_FASTQ | sed "s/.gz//")
        zcat $PREV_R1_FASTQ > $UNZIPPED_R1_FASTQ
        zcat $PREV_R1_FASTQ > $UNZIPPED_R2_FASTQ
        STAR --genomeDir \
            /oak/stanford/groups/cgawad/Reference_Files/GATK_Resource_Bundle_hg38/hg38_STAR_index/ \
            --runThreadN 2 --readFilesIn $UNZIPPED_R1_FASTQ $UNZIPPED_R2_FASTQ \
            --outFileNamePrefix $SAMPLE --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within --outSAMattributes Standard
        rm $UNZIPPED_R1_FASTQ $UNZIPPED_R2_FASTQ
        mv ${SAMPLE}Aligned.sortedByCoord.out.bam ${SAMPLE}_${REF_NAME}_aligned.bam
        echo "### Aligning RNA fastqs to $REF_NAME ### - END: $(date)"
    else
        echo "### Aligning DNA fastqs to $REF_NAME ### - START: $(date)"
        echo -e "Ref fasta: $REF_FASTA\nR1 fastq: $PREV_R1_FASTQ\nR2 fastq: $PREV_R2_FASTQ"
        bwa aln -t 2 $REF_FASTA $PREV_R1_FASTQ > ${SAMPLE}_${REF_NAME}_R1.sai
        bwa aln -t 2 $REF_FASTA $PREV_R2_FASTQ > ${SAMPLE}_${REF_NAME}_R2.sai
        bwa sampe -a 700 $REF_FASTA ${SAMPLE}_${REF_NAME}_R1.sai ${SAMPLE}_${REF_NAME}_R2.sai \
            $PREV_R1_FASTQ $PREV_R2_FASTQ | samtools view -b - | samtools sort -o ${SAMPLE}_${REF_NAME}_aligned.bam -
        samtools index ${SAMPLE}_${REF_NAME}_aligned.bam
        rm ${SAMPLE}_${REF_NAME}_R1.sai ${SAMPLE}_${REF_NAME}_R2.sai
        echo "### Aligning DNA fastqs to $REF_NAME ### - END: $(date)"
    fi
    
    echo "### Collecting $REF_NAME alignment metrics ### - START: $(date)"
    gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xmx32g" CollectAlignmentSummaryMetrics \
        -R $REF_FASTA -I ${SAMPLE}_${REF_NAME}_aligned.bam -O temp_${SAMPLE}_${REF_NAME}_ref_alignment_metrics.tsv
    echo "### Collecting $REF_NAME alignment metrics ### - END: $(date)"
    
    echo "### Filtering unmapped reads from $REF_NAME into a new BAM ### - START: $(date)"
    samtools view -b -f 4 ${SAMPLE}_${REF_NAME}_aligned.bam > ${SAMPLE}_no_${REF_NAME}.bam
    echo "### Filtering unmapped reads from $REF_NAME into a new BAM ### - END: $(date)"
    
    echo "### Converting unmapped reads from $REF_NAME from BAM to fastq ### - START: $(date)"
    gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xmx32g" SamToFastq -I ${SAMPLE}_no_${REF_NAME}.bam \
        -F ${SAMPLE}_no_${REF_NAME}${R1_SUFFIX} -F2 ${SAMPLE}_no_${REF_NAME}${R2_SUFFIX} --VALIDATION_STRINGENCY SILENT
    echo "### Converting unmapped reads from $REF_NAME from BAM to fastq ### - START: $(date)"
    
    PREV_R1_FASTQ=${SAMPLE}_no_${REF_NAME}${R1_SUFFIX}
    PREV_R2_FASTQ=${SAMPLE}_no_${REF_NAME}${R2_SUFFIX}
    rm ${SAMPLE}_${REF_NAME}_aligned.bam ${SAMPLE}_${REF_NAME}_aligned.bam.bai
done

for ((REF_INDEX = 0 ; REF_INDEX < ${#REF_FASTA_ARRAY[@]} ; REF_INDEX++)); do
    REF_FASTA=${REF_FASTA_ARRAY[$REF_INDEX]}
    REF_NAME=${REF_NAME_ARRAY[$REF_INDEX]}
    
    NEXT_INDEX=$((REF_INDEX+1))
    if [ $NEXT_INDEX -eq ${#REF_FASTA_ARRAY[@]} ]; then
        cp ${SAMPLE}_no_${REF_NAME}${R1_SUFFIX} ${SAMPLE}_ref_filtered${R1_SUFFIX}
        cp ${SAMPLE}_no_${REF_NAME}${R2_SUFFIX} ${SAMPLE}_ref_filtered${R2_SUFFIX}
        cp ${SAMPLE}_no_${REF_NAME}.bam ${SAMPLE}_ref_filtered.bam
        cp ${SAMPLE}_no_${REF_NAME}.bam.bai ${SAMPLE}_ref_filtered.bam.bai
    else
        rm ${SAMPLE}_no_${REF_NAME}${R1_SUFFIX} ${SAMPLE}_no_${REF_NAME}${R2_SUFFIX}
        rm ${SAMPLE}_no_${REF_NAME}.bam ${SAMPLE}_no_${REF_NAME}.bam.bai
    fi
done

if [ ! -f ${SAMPLE}_ref_filtered.bam ]; then
    echo "${SAMPLE}_ref_filtered.bam not found. Exiting with code 1"
    exit 1
fi
if [ $SKIP_TRIMMOMATIC -eq 0 ]; then
    rm ${SAMPLE}_trimmomatic_log.txt
    rm $R1_FASTQ $R2_FASTQ
    rm $UNPAIRED_R1_FASTQ $UNPAIRED_R2_FASTQ
fi
echo -e "END: $(date)\nRuntime: $(($(date +%s)-$START_TIME)) seconds"
