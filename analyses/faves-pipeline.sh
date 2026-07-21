tmp_dir=.

dm6_fasta="dm6.fasta"
dm6_prefix="dm6"
dm6_sim_prefix="dm6_prefix"
dm6_sim_reads="dm6.fastq.gz"
dm6_sim_gold="dm6.bed"

ecoliK12_fasta="ecoliK12.fasta"
ecoliK12_prefix="ecoliK12"
ecoliK12_sim_prefix="ecoliK12_prefix"
ecoliK12_sim_reads="ecoliK12.fastq.gz"
ecoliK12_sim_gold="ecoliK12.bed"

mtub_fasta="mtub.fasta"
mtub_prefix="mtub"
mtub_sim_prefix="mtub_prefix"
mtub_sim_reads="mtub.fastq.gz"
mtub_sim_gold="mtub.bed"

pf3D7_fasta="pf3D7.fasta"
pf3D7_prefix="pf3D7"
pf3D7_sim_prefix="pf3D7_prefix"
pf3D7_sim_reads="pf3D7.fastq.gz"
pf3D7_sim_gold="pf3D7.bed"

sacCer3_fasta="sacCer3.fasta"
sacCer3_prefix="sacCer3"
sacCer3_sim_prefix="sacCer3_prefix"
sacCer3_sim_reads="sacCer3.fastq.gz"
sacCer3_sim_gold="sacCer3.bed"

human_fasta="human.fasta"
human_prefix="human_prefix"
hg002_prefix="hg002_prefix"
hg002_reads="hg002.fastq.gz"
hg002_gold="HG002_GRCh38_1_22_v4.2.1_benchmark.bed"

datasets=(dm6 ecoliK12 mtub pf3D7 sacCer3 human)
tools=(fvs mm2-gtk) # e2i)

FAVES="../faves"

CONFIG_FILE="faves-config.tmp.sh"

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
# Portable replacement for `bedtools intersect -a A -b B [-wa -u | -v] | wc -l`.
# mode=u counts entries in A that overlap ANY entry in B (deduplicated by line).
# mode=v counts entries in A that overlap NO entry in B.
# Half-open intervals [start, end), same convention as BED.
bed_intersect_count () {
    local mode=$1
    local a=$2
    local b=$3
    awk -v mode="$mode" -F'\t' '
        FNR == NR {
            chr = $1
            n[chr] = (chr in n ? n[chr] : 0) + 1
            S[chr, n[chr]] = $2 + 0
            E[chr, n[chr]] = $3 + 0
            next
        }
        {
            astart = $2 + 0; aend = $3 + 0
            ov = 0; cnt = (($1 in n) ? n[$1] : 0)
            for (k = 1; k <= cnt; k++) {
                if (S[$1, k] < aend && E[$1, k] > astart) { ov = 1; break }
            }
            if ((mode == "u" && ov) || (mode == "v" && !ov)) c++
        }
        END { print c + 0 }
    ' "$b" "$a"
}

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
    local fasta=$1
    local ref_prefix=$2
    local sample_prefix=$3
    shift 3
    local fastq=("$@")

    if [ ! -f "$fasta" ]; then
        echo "Missing fasta $fasta"
        return
    fi

    for f in "${fastq[@]}"; do
        if [ ! -f "$f" ]; then
            echo "Missing reads $f"
            return
        fi
    done

    echo "Running FAVeS $ref_prefix ..."
    
    mkdir -p $tmp_dir/faves

    /bin/time -v ${FAVES} \
        -f $fasta \
        -q "${fastq[@]}" \
        -t 64 \
        -o $tmp_dir/faves/$sample_prefix.faves.snps.bed \
        -v 2> $tmp_dir/faves/$sample_prefix.faves.log
}

minimap2_pipeline() {
    local fasta=$1
    local ref_prefix=$2
    local sample_prefix=$3
    shift 3
    local fastq=("$@")

    echo "Running minimap2 $ref_prefix ..."

    if [ ! -f "$fasta" ]; then
        echo "Missing fasta $fasta"
        return
    fi

    for f in "${fastq[@]}"; do
        if [ ! -f "$f" ]; then
            echo "Missing reads $f"
            return
        fi
    done

    /bin/time -v minimap2 \
        -ax sr \
        -R '@RG\tID:1\tSM:hg002\tPL:ILLUMINA\tLB:lib1\tPU:unit1' \
        $fasta \
        "${fastq[@]}" \
        -t 64 > $tmp_dir/$sample_prefix.mm2.sam 2> $tmp_dir/faves/$sample_prefix.mm2.log

    /bin/time -v samtools sort $tmp_dir/$sample_prefix.mm2.sam > $tmp_dir/$sample_prefix.mm2.bam 2>> $tmp_dir/faves/$sample_prefix.mm2.log

    rm -f $tmp_dir/$sample_prefix.mm2.sam

    /bin/time -v samtools index $tmp_dir/$sample_prefix.mm2.bam 2>> $tmp_dir/faves/$sample_prefix.mm2.log
}

gatk_pipeline() {
    local fasta=$1
    local ref_prefix=$2
    local sample_prefix=$3
    local dict="${fasta%.fasta}.dict"

    if ! command -v gatk >/dev/null 2>&1; then
        echo "GATK not found, loading conda env..."
        conda activate gatk
        GATK_ENV_ACTIVATED=1
    fi

    echo "Running GATK $ref_prefix ..."

    if [ ! -f "$fasta" ]; then
        echo "Missing fasta $fasta"
        return
    fi

    if [ ! -f "$tmp_dir/$sample_prefix.mm2.bam" ]; then
        echo "Missing bam $tmp_dir/$sample_prefix.mm2.bam"
        return
    fi

    rm -f $tmp_dir/faves/$sample_prefix.gatk.log

    if [ ! -f  "$dict" ]; then 
        echo "Missing dict $dict . Creating ..." 
        gatk CreateSequenceDictionary -R $fasta -O $dict >> $tmp_dir/faves/$sample_prefix.gatk.log 2>> $tmp_dir/faves/$sample_prefix.gatk.log
    fi

    mkdir -p $tmp_dir/faves

    /bin/time -v gatk --java-options "-Xmx16g" HaplotypeCaller \
        --reference $fasta \
        --input $tmp_dir/$sample_prefix.mm2.bam \
        --output $tmp_dir/faves/$sample_prefix.mm2.gatk.vcf.gz \
        --native-pair-hmm-threads 64 >> "$tmp_dir/faves/$sample_prefix.gatk.log" 2>> "$tmp_dir/faves/$sample_prefix.gatk.log"
        
    /bin/time -v gatk SelectVariants \
        -R $fasta \
        -V $tmp_dir/faves/$sample_prefix.mm2.gatk.vcf.gz \
        --select-type-to-include SNP \
        -O $tmp_dir/faves/$sample_prefix.mm2.gatk.snps.vcf >> "$tmp_dir/faves/$sample_prefix.gatk.log" 2>> "$tmp_dir/faves/$sample_prefix.gatk.log"

    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' $tmp_dir/faves/$sample_prefix.mm2.gatk.snps.vcf | \
        awk 'BEGIN{OFS="\t"}{start=$2-1; end=start+length($4); print $1,start,end,$3,$4,$5}' > $tmp_dir/faves/$sample_prefix.mm2.gatk.snps.bed
}

ebwt2InDel_pipeline() {
    local fasta=$1
    local ref_prefix=$2
    local sample_prefix=$3
    shift 3
    local fastq=("$@")

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

for dataset in "${datasets[@]}"; do
    case "$dataset" in
        dm6)
            fasta=$dm6_fasta
            prefix=$dm6_prefix
            sim_prefix=$dm6_sim_prefix
            gold=$dm6_sim_gold
            reads=("${dm6_sim_reads[@]}")
            ;;
        ecoliK12)
            fasta=$ecoliK12_fasta
            prefix=$ecoliK12_prefix
            sim_prefix=$ecoliK12_sim_prefix
            gold=$ecoliK12_sim_gold
            reads=("${ecoliK12_sim_reads[@]}")
            ;;
        mtub)
            fasta=$mtub_fasta
            prefix=$mtub_prefix
            sim_prefix=$mtub_sim_prefix
            gold=$mtub_sim_gold
            reads=("${mtub_sim_reads[@]}")
            ;;
        pf3D7)
            fasta=$pf3D7_fasta
            prefix=$pf3D7_prefix
            sim_prefix=$pf3D7_sim_prefix
            gold=$pf3D7_sim_gold
            reads=("${pf3D7_sim_reads[@]}")
            ;;
        sacCer3)
            fasta=$sacCer3_fasta
            prefix=$sacCer3_prefix
            sim_prefix=$sacCer3_sim_prefix
            gold=$sacCer3_sim_gold
            reads=("${sacCer3_sim_reads[@]}")
            ;;
        human)
            fasta=$human_fasta
            prefix=$human_prefix
            sim_prefix=$hg002_prefix
            gold=$hg002_gold
            reads=("${hg002_reads[@]}")
            ;;    
    esac

    for tool in "${tools[@]}"; do
        case "$tool" in
            mm2-gtk) 
                minimap2_pipeline "$fasta" "$prefix" "$sim_prefix" "${reads[@]}"
                gatk_pipeline "$fasta" "$prefix" "$sim_prefix"
                snp_eval "MM2 + GATK vs Gold - $dataset" "$tmp_dir/faves/$sim_prefix.mm2.gatk.snps.bed" "$gold"
                ;;
            fvs)
                faves_pipeline "$fasta" "$prefix" "$sim_prefix" "${reads[@]}"
                snp_eval "FAVES vs Gold - $dataset" "$tmp_dir/faves/$sim_prefix.faves.snps.bed" "$gold"
                ;;
            e2i)
                ebwt2InDel_pipeline "$fasta" "$prefix" "$sim_prefix" "${reads[@]}"
                snp_eval "E2I vs Gold - $dataset" "$tmp_dir/faves/$sim_prefix.ebwt2InDel.5.snps.bed" "$gold"
                ;;
        esac
    done
done

# conda cleanup
if [ "$GATK_ENV_ACTIVATED" -eq 1 ]; then
    echo "Deactivating GATK conda env"
    conda deactivate
fi