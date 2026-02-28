: "${RUN_FVS:=true}"
: "${RUN_MM2:=true}"
: "${RUN_E2I:=true}"

tmp_dir="${tmp_dir:-}"
ref_dir="${ref_dir:-}"

hg38_prefix="${hg38_prefix:-}"
hg38_chr22_prefix="${hg38_chr22_prefix:-}"
hg002_prefix="${hg002_prefix:-}"
hg002_chr22_prefix="${hg002_chr22_prefix:-}"

CONFIG_FILE="faves-config.sh"

if [[ -n "$CONFIG_FILE" && -f "$CONFIG_FILE" ]]; then
    echo "Loading config from $CONFIG_FILE"
    source "$CONFIG_FILE"
fi


# -----------------------------------------
# -----------------------------------------
# SNP Callers
# -----------------------------------------
# -----------------------------------------
faves_pipeline() {
    local ref_prefix=$1
    local sample_prefix=$2
    local fastq=$3

    /bin/time -v faves \
        -f $ref_dir/$ref_prefix.fasta \
        -q $fastq \
        -t 64 \
        -o $tmp_dir/$sample_prefix.faves.bed \
        -p -v
}

minimap2_pipeline() {
    local ref_prefix=$1
    local sample_prefix=$2
    local fastq=$3

    /bin/time -v minimap2 \
        -ax sr \
        -R '@RG\tID:1\tSM:hg002\tPL:ILLUMINA\tLB:lib1\tPU:unit1' \
        $ref_dir/$ref_prefix.fasta \
        $fastq \
        -t 64 > $tmp_dir/$sample_prefix.mm2.sam

    /bin/time -v samtools sort $tmp_dir/$sample_prefix.mm2.sam > $tmp_dir/$sample_prefix.mm2.bam

    rm -f $tmp_dir/$sample_prefix.mm2.sam

    /bin/time -v samtools index $tmp_dir/$sample_prefix.mm2.bam

    if [ ! -f  "~/data/reference/$ref_prefix.dict" ]; then 
        gatk CreateSequenceDictionary -R $ref_dir/$ref_prefix.fasta -O $ref_dir/$ref_prefix.dict
    fi

    /bin/time -v gatk --java-options "-Xmx16g" HaplotypeCaller \
        --reference $ref_dir/$ref_prefix.fasta \
        --input $tmp_dir/$sample_prefix.mm2.bam \
        --output $tmp_dir/$sample_prefix.mm2.vcf.gz \
        --QUIET true \
        --native-pair-hmm-threads 64
    
    /bin/time -v gatk SelectVariants \
        -R $ref_dir/$ref_prefix.fasta \
        -V $tmp_dir/$sample_prefix.mm2.vcf.gz \
        --select-type-to-include SNP \
        -O $tmp_dir/$sample_prefix.mm2.snps.vcf

    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' $tmp_dir/$sample_prefix.mm2.snps.vcf | \
        awk 'BEGIN{OFS="\t"}{start=$2-1; end=start+length($4); print $1,start,end,$3,$4,$5}' > $tmp_dir/$sample_prefix.mm2.snps.bed
}

ebwt2InDel_pipeline() {
    local ref_prefix=$1
    local sample_prefix=$2
    local fasta=$3

    /bin/time -v BCR\_LCP\_GSA $fasta $tmp_dir/$sample_prefix.ebwt2InDel.bwt 1024

    /bin/time -v ebwt2InDel -1 $tmp_dir/$sample_prefix.ebwt2InDel.bwt -o $tmp_dir/$sample_prefix.ebwt2InDel.snp

    filter_snp $tmp_dir/$sample_prefix.ebwt2InDel.snp 5 >  $tmp_dir/$sample_prefix.ebwt2InDel.5.snp
}

source ~/.bashrc
conda activate gatk

# -----------------------------------------
# -----------------------------------------
# FAVeS
# -----------------------------------------
# -----------------------------------------
if [ "$RUN_FVS" = "true" ]; then
    # CHROMOSOME 22
    faves_pipeline $hg38_chr22_prefix $hg002_chr22_prefix ~/data/hg002/Illumina/HG002.GRCh38.chr22.2x250.fastq.gz

    # WGA
    faves_pipeline $hg38_prefix $hg002_prefix ~/data/hg002/Illumina/HG002.GRCh38.2x250.fastq.gz
fi

# -----------------------------------------
# -----------------------------------------
# minimap2
# -----------------------------------------
# -----------------------------------------
if [ "$RUN_MM2" = "true" ]; then
    # CHROMOSOME 22
    minimap2_pipeline $hg38_chr22_prefix $hg002_chr22_prefix ~/data/hg002/Illumina/HG002.GRCh38.chr22.2x250.fastq.gz

    # WGA
    minimap2_pipeline $hg38_prefix $hg002_prefix ~/data/hg002/Illumina/HG002.GRCh38.2x250.fastq.gz
fi

# -----------------------------------------
# -----------------------------------------
# ebwt2InDel
# -----------------------------------------
# -----------------------------------------
if [ "$RUN_E2I" = "true" ]; then
    # CHROMOSOME 22
    ebwt2InDel_pipeline $hg38_chr22_prefix $hg002_chr22_prefix ~/data/hg002/Illumina/HG002.GRCh38.chr22.2x250.fastq.gz

    # WGA
    ebwt2InDel_pipeline $hg38_prefix $hg002_prefix ~/data/hg002/Illumina/HG002.GRCh38.2x250.fastq.gz
fi

# conda cleanup
conda deactivate