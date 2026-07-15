#! /bin/bash

# please install simuG and make alias similar to:
# alias simuG="perl simuG.pl"

PROFILE="${PROFILE:-MSv3}"      # ART error model: MSv3 = MiSeq v3, supports up to 250 bp
READLEN="${READLEN:-250}"       # read length (bp)
FRAG_MEAN="${FRAG_MEAN:-500}"   # fragment/insert mean
FRAG_SD="${FRAG_SD:-50}"        # fragment/insert stddev
COVERAGE="${COVERAGE:-30}"      # simulated fold depth

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$ROOT/../data"

download() {
	local URL=$1
	local FILE_GZ=$2
	local FILE="${FILE_GZ%.gz}"

	wget $URL/$FILE_GZ
	gunzip $FILE_GZ
	mv $FILE "${FILE%.fna}.fasta"
	samtools faidx "${FILE%.fna}.fasta"
}

build_hap() {
	local FASTA=$1 PREFIX=$2 SNP_COUNT=$3 INDEL_COUNT=$4 CNV_COUNT=$5 INV_COUNT=$6 TRL_COUNT=$7;

	simuG \
		-refseq $FASTA \
		-snp_count "$SNP_COUNT" \
		-indel_count "$INDEL_COUNT" \
		-cnv_count "$CNV_COUNT" \
		-inversion_count "$INV_COUNT" \
		-translocation_count "$TRL_COUNT" \
		-titv_ratio 2.0 \
		-ins_del_ratio 1.0 \
		-indel_size_powerlaw_alpha 2.0 \
		-seed 100 \
		-prefix "$PREFIX"_sim;

	grep -v "^#" "${PREFIX}_sim.refseq2simseq.SNP.vcf" | awk 'BEGIN{OFS="\t"}{print $1, $2-1, $2, $4, $5}' > "${PREFIX}_sim.refseq2simseq.SNP.bed";
	grep -v "^#" "${PREFIX}_sim.refseq2simseq.INDEL.vcf" | awk 'BEGIN{OFS="\t"}{print $1, $2-1, $2, $4, $5}' > "${PREFIX}_sim.refseq2simseq.INDEL.bed";
}

build_dip() {
    local FASTA=$1 PREFIX=$2 SNP_COUNT=$3 INDEL_COUNT=$4 CNV_COUNT=$5 INV_COUNT=$6 TRL_COUNT=$7

    # draw THREE independent sets, ALL in original REF coords
    simuG -refseq "$FASTA" \
        -snp_count "$SNP_COUNT" -indel_count "$INDEL_COUNT" \
        -titv_ratio 2.0 -ins_del_ratio 1.0 -indel_size_powerlaw_alpha 2.0 \
        -seed 100 -prefix "${PREFIX}_shared"      # HOMOZYGOUS layer (both strands)

    simuG -refseq "$FASTA" \
        -snp_count "$SNP_COUNT" -indel_count "$INDEL_COUNT" \
        -titv_ratio 2.0 -ins_del_ratio 1.0 -indel_size_powerlaw_alpha 2.0 \
        -seed 201 -prefix "${PREFIX}_hap1priv"    # HET on hap1 only

    simuG -refseq "$FASTA" \
        -snp_count "$SNP_COUNT" -indel_count "$INDEL_COUNT" \
        -titv_ratio 2.0 -ins_del_ratio 1.0 -indel_size_powerlaw_alpha 2.0 \
        -seed 202 -prefix "${PREFIX}_hap2priv"    # HET on hap2 only

    # drop any private variant that collides with a shared
    # position, then merge shared+private and build each haplotype in ONE run
    drop_collisions() {
        awk 'FNR==NR{ if($0!~/^#/) s[$1"_"$2]=1; next }
             /^#/{ print; next }
             !(($1"_"$2) in s){ print }' "$1" "$2"
    }
    merge_vcf() {
        { grep "^#" "$1"; { grep -v "^#" "$1"; drop_collisions "$1" "$2" | grep -v "^#"; } | sort -k1,1 -k2,2n; } > "$3"
    }

    for T in SNP INDEL; do
        merge_vcf "${PREFIX}_shared.refseq2simseq.${T}.vcf" "${PREFIX}_hap1priv.refseq2simseq.${T}.vcf" "${PREFIX}_hap1.combined.${T}.vcf"
        merge_vcf "${PREFIX}_shared.refseq2simseq.${T}.vcf" "${PREFIX}_hap2priv.refseq2simseq.${T}.vcf" "${PREFIX}_hap2.combined.${T}.vcf"
    done

    simuG -refseq "$FASTA" -snp_vcf "${PREFIX}_hap1.combined.SNP.vcf" -indel_vcf "${PREFIX}_hap1.combined.INDEL.vcf" -prefix "${PREFIX}_hap1"
    simuG -refseq "$FASTA" -snp_vcf "${PREFIX}_hap2.combined.SNP.vcf" -indel_vcf "${PREFIX}_hap2.combined.INDEL.vcf" -prefix "${PREFIX}_hap2"

    # SINGLE truth BED per type, with genotype. All REF coords
    build_bed() {
        {
          grep -v "^#" "$1" | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$4,$5,"1/1"}'
          drop_collisions "$1" "$2" | grep -v "^#" | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$4,$5,"0/1"}'
          drop_collisions "$1" "$3" | grep -v "^#" | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$4,$5,"1/0"}'
        } | sort -k1,1 -k2,2n > "$4"
    }

    build_bed "${PREFIX}_shared.refseq2simseq.SNP.vcf"   "${PREFIX}_hap1priv.refseq2simseq.SNP.vcf"   "${PREFIX}_hap2priv.refseq2simseq.SNP.vcf"   "${PREFIX}_diploid.SNP.bed"
    build_bed "${PREFIX}_shared.refseq2simseq.INDEL.vcf" "${PREFIX}_hap1priv.refseq2simseq.INDEL.vcf" "${PREFIX}_hap2priv.refseq2simseq.INDEL.vcf" "${PREFIX}_diploid.INDEL.bed"
}

simulate_hap_reads() {
	local FASTA=$1 PREFIX=$2 COV=${3:-$COVERAGE} SEED=${4:-100};

	art_illumina -ss "$PROFILE" -i "$FASTA" -p -l "$READLEN" -f "$COV" -m "$FRAG_MEAN" -s "$FRAG_SD" -na -rs "$SEED" -o "${PREFIX}_R";
  	gzip -f "${PREFIX}_R1.fq" "${PREFIX}_R2.fq";
}

simulate_dip_reads() {
	local FASTA1=$1;
	local FASTA2=$2;
	local RATIO1=$3; # ex. 0.4
	local RATIO2=$4; # ex. 0.6

	local FASTA1=$1 FASTA2=$2 PREFIX=$3 RATIO1=$4 RATIO2=$5 COV=${6:-$COVERAGE}
	local COV1 COV2
	COV1=$(awk -v c="$COV" -v r="$RATIO1" 'BEGIN{printf "%.4f", c*r}')
	COV2=$(awk -v c="$COV" -v r="$RATIO2" 'BEGIN{printf "%.4f", c*r}')

	art_illumina -ss "$PROFILE" -i "$FASTA1" -p -l "$READLEN" -f "$COV1" -m "$FRAG_MEAN" -s "$FRAG_SD" -na -rs 201 -o "${PREFIX}_hap1_R";
	art_illumina -ss "$PROFILE" -i "$FASTA2" -p -l "$READLEN" -f "$COV2" -m "$FRAG_MEAN" -s "$FRAG_SD" -na -rs 202 -o "${PREFIX}_hap2_R";

 	awk '{if(NR%4==1)sub(/^@/,"@H1_"); print}' "${PREFIX}_hap1_R1.fq" >  "${PREFIX}_R1.fq";
	awk '{if(NR%4==1)sub(/^@/,"@H2_"); print}' "${PREFIX}_hap2_R1.fq" >> "${PREFIX}_R1.fq";
	awk '{if(NR%4==1)sub(/^@/,"@H1_"); print}' "${PREFIX}_hap1_R2.fq" >  "${PREFIX}_R2.fq";
	awk '{if(NR%4==1)sub(/^@/,"@H2_"); print}' "${PREFIX}_hap2_R2.fq" >> "${PREFIX}_R2.fq";
	gzip -f "${PREFIX}_R1.fq" "${PREFIX}_R2.fq";
	rm -f "${PREFIX}_hap1_R1.fq" "${PREFIX}_hap1_R2.fq" "${PREFIX}_hap2_R1.fq" "${PREFIX}_hap2_R2.fq";
}

mkdir -p "$ROOT/../data/dm6"  # diploid (2n=8)
cd "$ROOT/../data/dm6"
download "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT" "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz"
build_dip "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fasta" "dm6" 383000 48000 33 20 2
simulate_dip_reads "dm6_hap1.simseq.genome.fa" "dm6_hap2.simseq.genome.fa" "dm6" 0.4 0.6
mkdir simuG
mv dm6* simuG
cd "$ROOT"

mkdir -p "$ROOT/../data/ecoliK12"
cd "$ROOT/../data/ecoliK12"
download "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2" "GCF_000005845.2_ASM584v2_genomic.fna.gz"
build_hap "GCF_000005845.2_ASM584v2_genomic.fasta" "ecoliK12" 14000 1700 5 5 5
simulate_hap_reads "ecoliK12_sim.simseq.genome.fa" "ecoliK12"
mkdir simuG
mv ecoliK12* simuG
cd "$ROOT"

mkdir -p "$ROOT/../data/mtub"
cd "$ROOT/../data/mtub"
download "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2" "GCF_000195955.2_ASM19595v2_genomic.fna.gz"
build_hap "GCF_000195955.2_ASM19595v2_genomic.fasta" "mtub" 1300 160 3 3 3
simulate_hap_reads "mtub_sim.simseq.genome.fa" "mtub"
mkdir simuG
mv mtub* simuG
cd "$ROOT"

mkdir -p "$ROOT/../data/pf3D7"
cd "$ROOT/../data/pf3D7"
download "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/765/GCF_000002765.3_ASM276v1" "GCF_000002765.3_ASM276v1_genomic.fna.gz"
build_hap "GCF_000002765.3_ASM276v1_genomic.fasta" "pf3D7" 47000 5800 20 15 5
simulate_hap_reads "pf3D7_sim.simseq.genome.fa" "pf3D7"
mkdir simuG
mv pf3D7* simuG
cd "$ROOT"

mkdir -p "$ROOT/../data/sacCer3"
cd "$ROOT/../data/sacCer3"
download "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64" "GCF_000146045.2_R64_genomic.fna.gz"
build_hap "GCF_000146045.2_R64_genomic.fasta" "sacCer3" 60000 7500 20 15 5
simulate_hap_reads "sacCer3_sim.simseq.genome.fa" "sacCer3"
mkdir simuG
mv sacCer3* simuG
cd "$ROOT"

mkdir -p "$ROOT/../data/tair10"
cd "$ROOT/../data/tair10"   # diploid (2n=10)
download "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1" "GCF_000001735.4_TAIR10.1_genomic.fna.gz"
build_dip "GCF_000001735.4_TAIR10.1_genomic.fasta" "tair10" 225000 28000 27 17 2
simulate_dip_reads "tair10_hap1.simseq.genome.fa" "tair10_hap2.simseq.genome.fa" "tair10" 0.45 0.55
mkdir simuG
mv tair10* simuG
cd "$ROOT"

