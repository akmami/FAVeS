RUN_FVS_DM6="true"
RUN_MM2_DM6="true"
RUN_GTK_DM6="true"
RUN_E2I_DM6="true"
RUN_FVS_ECOLI="true"
RUN_MM2_ECOLI="true"
RUN_GTK_ECOLI="true"
RUN_E2I_ECOLI="true"
RUN_FVS_HUMAN="true"
RUN_MM2_HUMAN="true"
RUN_GTK_HUMAN="true"
RUN_E2I_HUMAN="true"

tmp_dir=.
ref_dir=.

dm6_prefix="dm6"
dm6_sim_prefix="dm6_prefix"
dm6_sim_reads="dm6.fastq.gz"
dm6_sim_gold="dm6.bed"

ecoli_prefix="ecoli"
ecoli_sim_prefix="ecoli_prefix"
ecoli_sim_reads="ecoli.fastq.gz"
ecoli_sim_gold="ecoli.bed"

hg38_chr22_prefix="hg38_chr22_prefix"
hg002_chr22_prefix="hg002_chr22_prefix"
hg002_chr22_reads="hg002.chr22.fastq.gz"
hg002_chr22_gold="HG002_GRCh38_22_v4.2.1_benchmark.bed"

hg38_prefix="hg38_prefix"
hg002_prefix="hg002_prefix"
hg002_reads="hg002.fastq.gz"
hg002_gold="HG002_GRCh38_1_22_v4.2.1_benchmark.bed"

CONFIG_FILE="faves-config.sh"

if [[ -n "$CONFIG_FILE" && -f "$CONFIG_FILE" ]]; then
    echo "Loading config from $CONFIG_FILE"
    source "$CONFIG_FILE"
fi

GATK_ENV_ACTIVATED=0


# -----------------------------------------
# -----------------------------------------
# Evaluation Callers
# -----------------------------------------
# -----------------------------------------
snp_eval () {
    local message=$1
    local calls=$2
    local gold=$3

    local TP=$(bedtools intersect -a "$calls" -b "$gold" -wa -u | wc -l)
    local FP=$(bedtools intersect -a "$calls" -b "$gold" -v | wc -l)
    local FN=$(bedtools intersect -a "$gold" -b "$calls" -v | wc -l)

    local PREC=$(echo "$TP / ($TP + $FP)" | bc -l)
    local REC=$(echo "$TP / ($TP + $FN)" | bc -l)
    local F1=$(echo "2 * $PREC * $REC / ($PREC + $REC)" | bc -l)

    echo $message
    echo "TP=$TP FP=$FP FN=$FN"
    echo "Precision=$PREC Recall=$REC F1=$F1"
    echo ""
}

# -----------------------------------------
# -----------------------------------------
# SNP Callers
# -----------------------------------------
# -----------------------------------------
faves_pipeline() {
    local ref_prefix=$1
    local sample_prefix=$2
    local fastq=$3

    echo "Running FAVeS $ref_prefix ..."
    
    mkdir -p $tmp_dir/faves

    /bin/time -v ../faves \
        -f $ref_dir/$ref_prefix.fasta \
        -q $fastq \
        -t 64 \
        -o $tmp_dir/faves/$sample_prefix.faves.snps.bed \
        -v 2> $tmp_dir/faves/$sample_prefix.faves.log
}

minimap2_pipeline() {
    local ref_prefix=$1
    local sample_prefix=$2
    local fastq=$3

    echo "Running minimap2 $ref_prefix ..."

    if [ ! -f "$ref_dir/$ref_prefix.fasta" ]; then
        echo "Missing fasta $ref_dir/$ref_prefix.fasta"
        return
    fi

    if [ ! -f "$fastq" ]; then
        echo "Missing reads $fastq"
        return
    fi

    /bin/time -v minimap2 \
        -ax sr \
        -R '@RG\tID:1\tSM:hg002\tPL:ILLUMINA\tLB:lib1\tPU:unit1' \
        $ref_dir/$ref_prefix.fasta \
        $fastq \
        -t 64 > $tmp_dir/$sample_prefix.mm2.sam 2> $tmp_dir/faves/$sample_prefix.mm2.log

    /bin/time -v samtools sort $tmp_dir/$sample_prefix.mm2.sam > $tmp_dir/$sample_prefix.mm2.bam 2>> $tmp_dir/faves/$sample_prefix.mm2.log

    rm -f $tmp_dir/$sample_prefix.mm2.sam

    /bin/time -v samtools index $tmp_dir/$sample_prefix.mm2.bam 2>> $tmp_dir/faves/$sample_prefix.mm2.log
}

gatk_pipeline() {
    local ref_prefix=$1
    local sample_prefix=$2

    if ! command -v gatk >/dev/null 2>&1; then
        echo "GATK not found, loading conda env..."
        source /home/akmuhammet/scripts/activate_conda
        conda activate gatk
        GATK_ENV_ACTIVATED=1
    fi

    echo "Running GATK $ref_prefix ..."

    if [ ! -f "$ref_dir/$ref_prefix.fasta" ]; then
        echo "Missing fasta $ref_dir/$ref_prefix.fasta"
        return
    fi

    if [ ! -f "$tmp_dir/$sample_prefix.mm2.bam" ]; then
        echo "Missing bam $tmp_dir/$sample_prefix.mm2.bam"
        return
    fi

    rm -f $tmp_dir/faves/$sample_prefix.gatk.log

    if [ ! -f  "~/data/reference/$ref_prefix.dict" ]; then 
        gatk CreateSequenceDictionary -R $ref_dir/$ref_prefix.fasta -O $ref_dir/$ref_prefix.dict 2>> $tmp_dir/faves/$sample_prefix.gatk.log
    fi

    mkdir -p $tmp_dir/faves

    /bin/time -v gatk --java-options "-Xmx16g" HaplotypeCaller \
        --reference $ref_dir/$ref_prefix.fasta \
        --input $tmp_dir/$sample_prefix.mm2.bam \
        --output $tmp_dir/faves/$sample_prefix.mm2.gatk.vcf.gz \
        --native-pair-hmm-threads 64 2>> $tmp_dir/faves/$sample_prefix.gatk.log
    
    /bin/time -v gatk SelectVariants \
        -R $ref_dir/$ref_prefix.fasta \
        -V $tmp_dir/faves/$sample_prefix.mm2.gatk.vcf.gz \
        --select-type-to-include SNP \
        -O $tmp_dir/faves/$sample_prefix.mm2.gatk.snps.vcf 2>> $tmp_dir/faves/$sample_prefix.gatk.log

    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' $tmp_dir/faves/$sample_prefix.mm2.gatk.snps.vcf | \
        awk 'BEGIN{OFS="\t"}{start=$2-1; end=start+length($4); print $1,start,end,$3,$4,$5}' > $tmp_dir/faves/$sample_prefix.mm2.gatk.snps.bed
    
    rm -f $ref_dir/$ref_prefix.dict
    rm -f $tmp_dir/faves/$sample_prefix.mm2.gatk.vcf.gz
    rm -f $tmp_dir/faves/$sample_prefix.mm2.gatk.snps.vcf
}

ebwt2InDel_pipeline() {
    local ref_prefix=$1
    local sample_prefix=$2
    local fasta=$3

    echo "Running ebwt2InDel $ref_prefix ..."

    if [ ! -f "$fasta" ]; then
        echo "Missing reads $fasta"
        return
    fi

    mkdir -p $tmp_dir/faves

    /bin/time -v BCR\_LCP\_GSA $fasta $tmp_dir/$sample_prefix.ebwt2InDel.bwt

    /bin/time -v ebwt2InDel -1 $tmp_dir/$sample_prefix.ebwt2InDel.bwt -o $tmp_dir/faves/$sample_prefix.ebwt2InDel.snp

    filter_snp $tmp_dir/faves/$sample_prefix.ebwt2InDel.snp 5 > $tmp_dir/faves/$sample_prefix.ebwt2InDel.5.snp
}

# -----------------------------------------
# -----------------------------------------
# FAVeS
# -----------------------------------------
# -----------------------------------------
if [ "$RUN_FVS_ECOLI" = "true" ]; then
    faves_pipeline $ecoli_prefix $ecoli_sim_prefix $ecoli_sim_reads
fi
if [ "$RUN_FVS_DM6" = "true" ]; then
    faves_pipeline $dm6_prefix $dm6_sim_prefix $dm6_sim_reads
fi
if [ "$RUN_FVS_HUMAN" = "true" ]; then
    faves_pipeline $hg38_chr22_prefix $hg002_chr22_prefix $hg002_chr22_reads
    faves_pipeline $hg38_prefix $hg002_prefix $hg002_reads
fi

# Evaluation
if [ "$RUN_FVS_ECOLI" = "true" ]; then
    snp_eval "FAVES vs Gold - ecoli" "$tmp_dir/faves/$ecoli_sim_prefix.faves.snps.bed" "$ecoli_sim_gold"
fi
if [ "$RUN_FVS_DM6" = "true" ]; then
    snp_eval "FAVES vs Gold - dm6" "$tmp_dir/faves/$dm6_sim_prefix.faves.snps.bed" "$dm6_sim_gold"
fi
if [ "$RUN_FVS_HUMAN" = "true" ]; then
    snp_eval "FAVES vs Gold - hg002.chr22" "$tmp_dir/faves/$hg002_chr22_prefix.faves.snps.bed" "$hg002_chr22_gold"
    snp_eval "FAVES vs Gold - hg002" "$tmp_dir/faves/$hg002_prefix.faves.snps.bed" "$hg002_gold"
fi

# -----------------------------------------
# -----------------------------------------
# minimap2
# -----------------------------------------
# -----------------------------------------
if [ "$RUN_MM2_ECOLI" = "true" ]; then
    minimap2_pipeline $ecoli_prefix $ecoli_sim_prefix $ecoli_sim_reads
fi
if [ "$RUN_MM2_DM6" = "true" ]; then
    minimap2_pipeline $dm6_prefix $dm6_sim_prefix $dm6_sim_reads
fi
if [ "$RUN_MM2_HUMAN" = "true" ]; then
    minimap2_pipeline $hg38_chr22_prefix $hg002_chr22_prefix $hg002_chr22_reads
    minimap2_pipeline $hg38_prefix $hg002_prefix $hg002_reads
fi

# -----------------------------------------
# -----------------------------------------
# GATK
# -----------------------------------------
# -----------------------------------------
if [ "$RUN_GTK_ECOLI" = "true" ]; then
    gatk_pipeline $ecoli_prefix $ecoli_sim_prefix
fi
if [ "$RUN_GTK_DM6" = "true" ]; then
    gatk_pipeline $dm6_prefix $dm6_sim_prefix
fi
if [ "$RUN_GTK_HUMAN" = "true" ]; then
    gatk_pipeline $hg38_chr22_prefix $hg002_chr22_prefix
    gatk_pipeline $hg38_prefix $hg002_prefix
fi

if [ "$RUN_GTK_ECOLI" = "true" ]; then
    snp_eval "MM2 + GATK vs Gold - ecoli" "$tmp_dir/faves/$ecoli_sim_prefix.mm2.gatk.snps.bed" "$ecoli_sim_gold"
fi
if [ "$RUN_GTK_DM6" = "true" ]; then
    snp_eval "MM2 + GATK vs Gold - dm6" "$tmp_dir/faves/$dm6_sim_prefix.mm2.gatk.snps.bed" "$dm6_sim_gold"
fi
if [ "$RUN_GTK_HUMAN" = "true" ]; then
    snp_eval "MM2 + GATK vs Gold - hg002.chr22" "$tmp_dir/faves/$hg002_chr22_prefix.mm2.gatk.snps.bed" "$hg002_chr22_gold"
    snp_eval "MM2 + GATK vs Gold - hg002" "$tmp_dir/faves/$hg002_prefix.mm2.gatk.snps.bed" "$hg002_gold"
fi

# -----------------------------------------
# -----------------------------------------
# ebwt2InDel
# -----------------------------------------
# -----------------------------------------
if [ "$RUN_E2I_ECOLI" = "true" ]; then
    ebwt2InDel_pipeline $ecoli_prefix $ecoli_sim_prefix $ecoli_sim_reads
fi
if [ "$RUN_E2I_DM6" = "true" ]; then
    ebwt2InDel_pipeline $dm6_prefix $dm6_sim_prefix $dm6_sim_reads
fi
if [ "$RUN_E2I_HUMAN" = "true" ]; then
    ebwt2InDel_pipeline $hg38_chr22_prefix $hg002_chr22_prefix $hg002_chr22_reads
    ebwt2InDel_pipeline $hg38_prefix $hg002_prefix $hg002_reads
fi

# Evaluation
if [ "$RUN_E2I_ECOLI" = "true" ]; then
    snp_eval "E2I vs Gold - ecoli" "$tmp_dir/faves/$ecoli_sim_prefix.ebwt2InDel.5.snps.bed" "$ecoli_sim_gold"
fi
if [ "$RUN_E2I_DM6" = "true" ]; then
    snp_eval "E2I vs Gold - dm6" "$tmp_dir/faves/$dm6_sim_prefix.ebwt2InDel.5.snps.bed" "$dm6_sim_gold"
fi
if [ "$RUN_E2I_HUMAN" = "true" ]; then
    snp_eval "E2I vs Gold - hg002.chr22" "$tmp_dir/faves/$hg002_chr22_prefix.ebwt2InDel.5.snps.bed" "$hg002_chr22_gold"
    snp_eval "E2I vs Gold - hg002" "$tmp_dir/faves/$hg002_prefix.ebwt2InDel.5.snps.bed" "$hg002_gold"
fi

# conda cleanup
if [ "$GATK_ENV_ACTIVATED" -eq 1 ]; then
    echo "Deactivating GATK conda env"
    conda deactivate
fi