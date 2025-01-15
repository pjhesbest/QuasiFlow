reference=/home/phesketh/Documents/GitHub/QuasiFlow/db/K03455.1_HIV/K03455.1_HIV.fasta
reference_index=/home/phesketh/Documents/GitHub/QuasiFlow/db/K03455.1_HIV/K03455.1
nthreads=2
sampleID=Prueba2-6

# Make the BWA-MEM2 index
    bwa-mem2 index ${reference} -p Chromosome

# Map to BWA-MEM2 database
    bwa-mem2 mem -t 2 ${reference_index} \
                ${sampleID}/${sampleID}_val_1.fq.gz \
                ${sampleID}/${sampleID}_val_2.fq.gz \
                > ${sampleID}/${sampleID}.sam

# sort .sam file and convert to .bam file
    samtools view -bS ${sampleID}/${sampleID}.sam \
                | samtools sort - -o ${sampleID}/${sampleID}.bam

    rm ${sampleID}/${sampleID}.sam

# call variants
    bcftools mpileup -Ou -f ${reference} ${sampleID}/${sampleID}.bam | bcftools call -mv -Oz -o ${sampleID}/${sampleID}.calls.vcf.gz

    bcftools index ${sampleID}/${sampleID}.calls.vcf.gz

# normalize indels
    bcftools norm -f ${reference} \
            ${sampleID}/${sampleID}.calls.vcf.gz \
            -Ob -o ${sampleID}/${sampleID}.calls.norm.bcf

# filter adjacent indels within 5bp
    bcftools filter --IndelGap 5 \
            ${sampleID}/${sampleID}.calls.norm.bcf \
            -Ob -o ${sampleID}/${sampleID}.calls.norm.flt-indels.bcf

# apply variants to create consensus sequence
    cat ${reference} \
            | bcftools consensus ${sampleID}/${sampleID}.calls.vcf.gz \
            > ${sampleID}/${sampleID}.consensus.fa

# output IUPAC ambiguity codes based on REF+ALT columns (regardless of genotype)
    cat ${reference} \
            | bcftools consensus \
            --iupac-codes ${sampleID}/${sampleID}.calls.vcf.gz \
            > ${sampleID}/${sampleID}.consensus.fa

# output IUPAC ambiguity codes based on sample genotypes
    cat ${reference} | bcftools consensus \
            --haplotype I ${sampleID}/${sampleID}.calls.vcf.gz \
            > ${sampleID}/${sampleID}.consensus.fa


for sampleID in Prueba*; do

> samtools sort ${sampleID}/${sampleID}.sam > ${sampleID}/${sampleID}.bam
> rm ${sampleID}/${sampleID}.sam
> samtools index ${sampleID}/${sampleID}.bam
> done


# cannot be performed on gzipped files
for sampleID in Prueba*; do

cd ${sampleID}

    gunzip ${sampleID}_val_1.fq.gz; gunzip ${sampleID}_val_2.fq.gz

quasitools hydra ${sampleID}_val_1.fq \
                ${sampleID}_val_2.fq -o . \
                --generate_consensus \
                --reporting_threshold 1 \
                --consensus_pct 20 \
                --length_cutoff 250 \
                --score_cutoff 20 \
                --min_variant_qual 20 \
                --min_dp 100 \
                --min_freq 0.01 \
                --min_ac 5

virulign HIV-HXB2-pol.xml HIVdb.fasta \
            --exportKind  GlobalAlignment \
            --exportAlphabet Nucleotides \
            --exportWithInsertions no \
            --exportReferenceSequence yes \
            --nt-debug Failed > HIVrt.fasta

    pigz --best ${sampleID}_val_1.fq ${sampleID}_val_2.fq

sierrapy --virus HIV1 fasta consensus.fasta -o ${sampleID}.sierrapy.hiv1.json

python3 /home/phesketh/Documents/projects/HIV/data/sierra-local/scripts/json2csv.py ${sampleID}.sierrapy.hiv1*json ${sampleID}.sierrapy.hiv1.csv

cd ../

done


