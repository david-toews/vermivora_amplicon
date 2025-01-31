
##Build reference amplicon scaffold
bowtie2-build /go/to/folder/GWBW_amplicon_test/amplicon_reference_4a_20_24_25_z1_z2.fa /go/to/folder/GWBW_amplicon_test/amplicon_reference_4a_20_24_25_z1_z2

java -jar /storage/group/dut374/default/sw/picard-3.1.0/bin/picard.jar CreateSequenceDictionary R=/go/to/folder/GWBW_amplicon_test/amplicon_reference_4a_20_24_25_z1_z2.fa O=/go/to/folder/GWBW_amplicon_test/amplicon_reference_4a_20_24_25_z1_z2.dict

samtools faidx /go/to/folder/GWBW_amplicon_test/amplicon_reference_4a_20_24_25_z1_z2.fa


# Navigate to the working directory
cd /go/to/folder

# Find all R1 fastq.gz files in the GWBW_amplicon_test directory
readpair_array=( $(find GWBW_amplicon_test -name "*R1_001.fastq.gz") )

# Loop through the found files
for j in "${readpair_array[@]}"; do
    # Remove "_L001_R1_001.fastq.gz" from the filename to get the sample name
    i=$(echo "$j" | sed "s/_L001_R1_001.fastq.gz//")
    sample=$(basename "$i")

    echo "Processing sample: $sample"

    # Create the SLURM job script
    cat <<EOF > "${sample}.align.sh"
#!/bin/bash
#SBATCH --account=######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=10:00:00
#SBATCH --mem=2GB
#SBATCH -o gwbw_amplicon_test.log
#SBATCH --job-name=gwbw_amplicon_test

### Load modules
module use /storage/group/dut374/default/sw/modules
module load all

nohup bowtie2 -p 1 --very-sensitive-local --local -N 0 --phred33 \
-x /go/to/folder/GWBW_amplicon_test/amplicon_reference_4a_20_24_25_z1_z2 \
--rg-id $sample \
--rg SM:$sample \
-1 /go/to/folder/GWBW_amplicon_test/${sample}_L001_R1_001.fastq.gz \
-2 /go/to/folder/GWBW_amplicon_test/${sample}_L001_R2_001.fastq.gz \
-X 100 \
-S /go/to/folder/GWBW_amplicon_test/BAM/$sample.sam
EOF
    # Submit the generated script to SLURM
    sbatch "./${sample}.align.sh"
done


###Compress, sort and index BAM files

# Set the main directory
cd /go/to/folder

# Find all SAM files in the GWBW_amplicon_test/BAM directory
sam_files=( $(find GWBW_amplicon_test/BAM -name "*.sam") )

# Loop through each SAM file
for sam_file in "${sam_files[@]}"; do
    # Extract the base name without extension
    sample=$(basename "$sam_file" .sam)

    echo "Processing sample: $sample"

    # Create the SLURM job script
    cat <<EOF > "${sample}.samtools.sh"

#!/bin/bash
#SBATCH --account=######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=4GB
#SBATCH -o ${sample}_samtools.log
#SBATCH --job-name=${sample}_samtools

### Load modules
module use /storage/group/dut374/default/sw/modules
module load samtools

### Step 1: Convert SAM to BAM
samtools view -bS /go/to/folder/GWBW_amplicon_test/BAM/${sample}.sam > /go/to/folder/GWBW_amplicon_test/BAM/${sample}.bam

### Step 2: Sort BAM
samtools sort /go/to/folder/GWBW_amplicon_test/BAM/${sample}.bam -o /go/to/folder/GWBW_amplicon_test/BAM/${sample}_sorted.bam

### Step 3: Index Sorted BAM
samtools index /go/to/folder/GWBW_amplicon_test/BAM/${sample}_sorted.bam
EOF

    # Submit the generated script to SLURM
    sbatch "./${sample}.samtools.sh"
done

###Run next portion in "interactive" mode
salloc --account=###### -N 1 -n 8 --mem-per-cpu=5000 -t 3:00:00


##Call SNPs with GATK

##make amplicon_sites.bed
cat <<EOF > amplicon_sites.bed
amplicon_reference_4a_20_24_25_z1_z2 227 228
amplicon_reference_4a_20_24_25_z1_z2 610 611
amplicon_reference_4a_20_24_25_z1_z2 1084 1085
amplicon_reference_4a_20_24_25_z1_z2 1554 1555
amplicon_reference_4a_20_24_25_z1_z2 1801 1802
amplicon_reference_4a_20_24_25_z1_z2 2263 2264
EOF


##Call SNPs with GATK: using -L to force it to call SNPs for just those sites defined above

# Define directories and reference
BAM_DIR="/go/to/folder/GWBW_amplicon_test/BAM"
VCF_DIR="/go/to/folder/GWBW_amplicon_test/VCF"
REFERENCE="/go/to/folder/GWBW_amplicon_test/amplicon_reference_4a_20_24_25_z1_z2.fa"

# Iterate over all sorted BAM files in the BAM directory
for BAM_FILE in ${BAM_DIR}/*sorted.bam; do
    # Get the base name of the BAM file (without path and extension)
    SAMPLE_NAME=$(basename ${BAM_FILE} .sorted.bam)

    # Output VCF file path
    VCF_OUTPUT="${VCF_DIR}/${SAMPLE_NAME}.g.vcf.gz"

    # Run HaplotypeCaller for each BAM file
    /storage/home/dut374/gatk-4.6.1.0/gatk HaplotypeCaller \
        -R ${REFERENCE} \
        -I ${BAM_FILE} \
        -O ${VCF_OUTPUT} \
        -ERC GVCF \
--ploidy 2 \
--max-alternate-alleles 2 \
--output-mode EMIT_ALL_CONFIDENT_SITES \
        -L /go/to/folder/amplicon_sites.bed
done



##The obvious next step here would be to combine ths GVCFs like normal and do joint genotyping, but GATK was doing weird
##stuff even with the combining of the g.vcfs. For instance, in some files where the individual genotypes looked fine
##the combined VCF coded them as missing. So in the end we just worked with the code below to extract
#the individual genotypes. Note-—they aren't polarized.

INPUT_DIR="/go/to/folder/GWBW_amplicon_test/VCF/"
OUTPUT_FILE="raw_genotypes.txt"

# Clear the output file if it exists and create the header
echo -e "CHROM\tPOS" > "$OUTPUT_FILE"

# Create a temporary header file
HEADER_FILE="header_temp.txt"
echo -ne "CHROM\tPOS" > "$HEADER_FILE"

# Process each gVCF file
for FILE in ${INPUT_DIR}/*.g.vcf.gz; do
    SAMPLE_NAME=$(basename "$FILE" .g.vcf.gz)  # Extract sample name
    echo "Processing $SAMPLE_NAME..."

    # Extract genotypes
    bcftools query -f '%CHROM\t%POS\t[%GT]\n' "$FILE" > "${SAMPLE_NAME}_genotypes.txt"

    # Store sample names as headers
    echo -ne "\t$SAMPLE_NAME" >> "$HEADER_FILE"

    # Merge results into the main output file
    if [[ $(wc -l < "$OUTPUT_FILE") -eq 1 ]]; then
        cp "${SAMPLE_NAME}_genotypes.txt" "$OUTPUT_FILE"
    else
        paste "$OUTPUT_FILE" <(cut -f3 "${SAMPLE_NAME}_genotypes.txt") > temp && mv temp "$OUTPUT_FILE"
    fi

    # Clean up
    rm "${SAMPLE_NAME}_genotypes.txt"
done

# Add header row to the final output
cat "$HEADER_FILE" > temp && echo "" >> temp && cat "$OUTPUT_FILE" >> temp && mv temp "$OUTPUT_FILE"
rm "$HEADER_FILE"