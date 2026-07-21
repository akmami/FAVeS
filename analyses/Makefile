SHELL := /bin/bash

# directories
CUR    := $(shell pwd)
REPO   := $(abspath $(CUR)/..)
SKETCH := $(REPO)/sketch
LIB    := $(REPO)/lib
BIN    := $(REPO)/bin
OUT    := $(REPO)/out
BUILD  := $(OUT)/.build

# toolchain
GXX      := g++
GCC      := gcc
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -Wpedantic
SIMD     := -mavx2 -msse4.1
# optional: use /usr/bin/time -v when available, otherwise run untimed
TIME     ?= $(shell command -v /usr/bin/time >/dev/null 2>&1 && echo '/usr/bin/time -v')
JOBS     ?= $(shell nproc 2>/dev/null || echo 4)
PYTHON   ?= python3
JOBS 	 ?= 8

# shortnames
dm6_NAME      = D mel.
ecoliK12_NAME = E. coli
mtub_NAME     = Mtb
pf3D7_NAME    = Pf
sacCer3_NAME  = Yeast
human_NAME    = Human

# genomes: <name>_FASTA. single: GENOME=<name>; batch: GENOMES="a b c"
dm6_FASTA      ?= ../data/dm6/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fasta
ecoliK12_FASTA ?= ../data/ecoliK12/GCF_000005845.2_ASM584v2_genomic.fasta
mtub_FASTA     ?= ../data/mtub/GCF_000195955.2_ASM19595v2_genomic.fasta
pf3D7_FASTA    ?= ../data/pf3D7/GCF_000002765.3_ASM276v1_genomic.fasta
sacCer3_FASTA  ?= ../data/sacCer3/GCF_000146045.2_R64_genomic.fasta
human_FASTA    ?= ../data/human/human_v38.fasta
GENOME  ?= human
GENOMES ?= ecoliK12 mtub sacCer3 pf3D7 dm6 human
FASTA   ?= $($(GENOME)_FASTA)
NAME    ?= $($(GENOME)_NAME)

# analysis knobs
KMAX      ?= 512
LHIST     ?= 2048
KSWEEP    ?= 11 13 15 17 19 21 23 25 27

# contiguous-k sketchers share the uniqueness-vs-k axis; strobemer (span) and lcp (level) use different parameters.
SWEEP_METHODS ?= kmer minimizer fracmh syncmer blend
SWEEP_GENOMES ?= $(GENOMES)
SWEEP_W       ?= 11      # minimizer window
SWEEP_MOD     ?= 6       # FracMinHash modulus
SWEEP_S       ?= 10      # syncmer s-mer size
SWEEP_B       ?= 50      # blend bits
SWEEP_N       ?= 5       # blend n

# method registry: add a technique = add one block
# <m>_MACRO  : -D define selecting the branch in sketch-seeds.cpp
# <m>_LIB    : static library base name (lib<LIB>.a)
# <m>_SRC    : source file relative to sketch/
# <m>_CC     : compiler for the library object (g++ for .cpp / kmer / lcp)
# <m>_PARAMS : run grid; tuples separated by '|', each is the args after FASTA
METHODS := kmer fracmh minimizer syncmer blend strobemer strobemer2 lcp

kmer_MACRO  := KMER
kmer_LIB    := kmer
kmer_SRC    := kmer.c
kmer_CC     := $(GXX)
kmer_PARAMS ?= 15 | 21

fracmh_MACRO  := FRACMH
fracmh_LIB    := fracmh
fracmh_SRC    := fracmh.c
fracmh_CC     := $(GCC)
fracmh_PARAMS ?= 15 $(SWEEP_MOD) | 21 $(SWEEP_MOD)

minimizer_MACRO  := MINIMIZER
minimizer_LIB    := minimizer
minimizer_SRC    := minimizer.c
minimizer_CC     := $(GCC)
minimizer_PARAMS ?= 21 11 | 15 10

syncmer_MACRO  := SYNCMER
syncmer_LIB    := syncmer
syncmer_SRC    := syncmer.c
syncmer_CC     := $(GCC)
syncmer_PARAMS ?= 21 10 | 15 10 | 15 5

blend_MACRO  := BLEND
blend_LIB    := blend
blend_SRC    := blend.c
blend_CC     := $(GCC)
blend_PARAMS ?= 21 11 50 3 | 21 11 50 5 | 21 11 50 7 | 15 10 50 3 | 15 10 50 5 | 15 10 50 7

strobemer_MACRO  := STROBEMER
strobemer_LIB    := strobemer
strobemer_SRC    := strobemer/strobemer.cpp
strobemer_CC     := $(GXX)
strobemer_PARAMS ?= 2 22 6 16

strobemer2_MACRO  := STROBEMER2
strobemer2_LIB    := strobemer2
strobemer2_SRC    := strobemer/strobemer2.cpp
strobemer2_CC     := $(GXX)
strobemer2_PARAMS ?= 250 17

lcp_MACRO  := LCP
lcp_LIB    := lcp
lcp_SRC    := lps.c
lcp_CC     := $(GXX)
lcp_PARAMS ?= 2 1 | 3 1 | 4 1

# libdivsufsort
DIVSUFSORT_PREFIX  ?= $(REPO)/third_party/libdivsufsort
DIVSUFSORT_URL     ?= https://github.com/y-256/libdivsufsort.git
DIVSUFSORT_TARBALL ?= https://github.com/y-256/libdivsufsort/archive/refs/heads/master.tar.gz
DS_INC  := $(if $(wildcard $(DIVSUFSORT_PREFIX)/include/divsufsort.h),-I$(DIVSUFSORT_PREFIX)/include,)
DS_LDIR := $(if $(wildcard $(DIVSUFSORT_PREFIX)/lib/libdivsufsort.a),-L$(DIVSUFSORT_PREFIX)/lib,)
DS_LIBS := $(DS_LDIR) -ldivsufsort -ldivsufsort64

.DEFAULT_GOAL := help
.PHONY: help dirs clean divsufsort matchlen matchlenAll sweep plots $(METHODS)
.SECONDARY:

# help
help:
	@echo "FAVeS analyses (single Makefile). Targets:"
	@echo "  make divsufsort                     build libdivsufsort locally"
	@echo "  make matchlen       [GENOME=human]  exact p(k)=CDF(L)            [1,2,3,5,7]"
	@echo "  make matchlenAll    [GENOMES=...]   match-length per genome      [4]"
	@echo "  make sketch-methods [GENOME=human]  run one technique's grid"
	@echo "  make sketch-all                     run all technique's grid"
	@echo "  make <m>                            shorthand for sketch-<m>"
	@echo "  make sweep          [GENOME=human]  kmer/minimizer/fracmh sweep  [6]"
	@echo "  make plots                          fit (alpha,beta) + figures"
	@echo "  methods : $(METHODS)"
	@echo "  genomes : $(GENOMES)"

dirs:
	@mkdir -p $(BIN) $(OUT) $(LIB) $(BUILD)

clean:
	@rm -rf $(BUILD) $(BIN)/sketch-seeds.* $(BIN)/match-length; echo "cleaned"

# libdivsufsort
divsufsort:
	@if [ -f "$(DIVSUFSORT_PREFIX)/lib/libdivsufsort64.a" ]; then \
		echo "libdivsufsort already present at $(DIVSUFSORT_PREFIX)"; \
		exit 0; \
	fi; \
	set -e; command -v cmake >/dev/null || { echo "ERROR: need cmake (try: module load cmake)"; exit 1; }; \
	w=$$(mktemp -d); \
	echo ">> fetching libdivsufsort"; \
	( git clone --depth 1 $(DIVSUFSORT_URL) $$w/src >/dev/null 2>&1 ) || \
		( curl -fsSL $(DIVSUFSORT_TARBALL) -o $$w/ds.tgz && mkdir -p $$w/src && tar xzf $$w/ds.tgz -C $$w/src --strip-components=1 ); \
	echo ">> building (static, 32+64-bit) -> $(DIVSUFSORT_PREFIX)"; \
	mkdir -p $$w/build; cd $$w/build; \
	cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_DIVSUFSORT64=ON -DBUILD_SHARED_LIBS=OFF \
		-DCMAKE_INSTALL_PREFIX=$(DIVSUFSORT_PREFIX) $$w/src >/dev/null; \
	make -j$(JOBS) >/dev/null; make install >/dev/null; rm -rf $$w; \
	echo ">> libdivsufsort installed"

# match-length
$(BIN)/match-length: match-length.cpp | dirs
	$(GXX) $(CXXFLAGS) $(DS_INC) -o $@ $< $(DS_LIBS)

matchlen: $(BIN)/match-length
	@mkdir -p $(OUT)/match-length
	$(TIME) $(BIN)/match-length "$(FASTA)" "$(NAME)" $(KMAX) $(LHIST) > $(OUT)/match-length/$(GENOME).json
	@echo "-> $(OUT)/match-length/$(GENOME).json"

matchlenAll: $(BIN)/match-length
	@mkdir -p $(OUT)/match-length; \
	declare -A F=( $(foreach g,$(GENOMES),[$(g)]="$($(g)_FASTA)") ); \
	declare -A S=( $(foreach g,$(GENOMES),[$(g)]="$($(g)_NAME)") ); \
	for g in $(GENOMES); do f=$${F[$$g]}; s=$${S[$$g]}; \
		if [ -f "$$f" ]; then \
			echo ">> match-length $$g ($$s)"; \
			$(BIN)/match-length "$$f" "$$s" $(KMAX) $(LHIST) > $(OUT)/match-length/$$g.json; \
		else \
			echo "  skip $$g ($$f not found)"; \
		fi; \
	done; \
	echo "-> $(OUT)/match-length/"

# generic sketch binary: builds lib + $(BIN)/sketch-seeds.<method>
$(BIN)/sketch-seeds.%: sketch-seeds.cpp | dirs
	@m=$*; macro=$($*_MACRO); lib=$($*_LIB); src=$($*_SRC); cc=$($*_CC); \
	if [ -z "$$macro" ]; then \
		echo "ERROR: unknown method '$$m' (know: $(METHODS))"; \
		exit 1; \
	fi; \
	std=""; case "$$src" in *.cpp) std="-std=c++17";; esac; \
	echo ">> lib$$lib.a  ($$src)"; \
	$$cc -O3 $(SIMD) $$std -c $(SKETCH)/$$src -o $(BUILD)/$$lib.o; \
	ar rcs $(LIB)/lib$$lib.a $(BUILD)/$$lib.o; \
	echo ">> sketch-seeds.$$m  (-D$$macro)"; \
	$(GXX) $(CXXFLAGS) $(SIMD) -D$$macro -c sketch-seeds.cpp -o $(BUILD)/ss.$$m.o; \
	$(GXX) $(CXXFLAGS) -o $@ $(BUILD)/ss.$$m.o -L$(LIB) -l$$lib

# run one technique over its parameter grid on $(GENOME)
sketch-%: $(BIN)/sketch-seeds.%
	@m=$*; d=$(OUT)/$$m; out=$$d/$(GENOME).json; \
	prog="$(BIN)/sketch-seeds.$$m \"$(FASTA)\" \"$(NAME)\""; \
	mkdir -p $$d; \
	echo ">> $$m x $(GENOME)"; \
	rm -f $$out $$d/tmp.$(GENOME).* $$d/$(GENOME).err; \
	i=0; \
	echo '$($*_PARAMS)' | tr '|' '\n' | sed '/^$$/d' | \
		while read -r args; do i=$$((i+1)); \
			echo "$$prog $$args > $$d/tmp.$(GENOME).$$i 2>> $$d/$(GENOME).err"; \
		done | \
		parallel -j $(JOBS); \
	{ echo "["; cat $$(ls -v $$d/tmp.$(GENOME).*); } > $$out; \
	perl -0777 -i -pe 's/,\s*\z/\n/' $$out; \
	echo "]" >> $$out; \
	rm -f $$d/tmp.$(GENOME).*; \
	[ -f "$$d/$(GENOME).err" ] && [ ! -s "$$d/$(GENOME).err" ] && rm "$$d/$(GENOME).err"; \
	echo "-> $$out"

# shorthand: `make kmer` == `make sketch-kmer`
$(METHODS): %: sketch-%

sketch-methods: $(addprefix sketch-,$(METHODS))
	@echo ">> all methods done -> $(OUT)"

sketch-all:
	@for g in $(GENOMES); do \
		echo "==================  $$g  =================="; \
		$(MAKE) -f $(firstword $(MAKEFILE_LIST)) --no-print-directory sketch-methods GENOME=$$g || exit $$?; \
	done

# sketching-invariance sweep
sweep: $(foreach m,$(SWEEP_METHODS),$(BIN)/sketch-seeds.$(m))
	@mkdir -p $(OUT)/sweep; \
	declare -A F=( $(foreach g,$(SWEEP_GENOMES),[$(g)]="$($(g)_FASTA)") ); \
	declare -A S=( $(foreach g,$(SWEEP_GENOMES),[$(g)]="$($(g)_NAME)") ); \
	for g in $(SWEEP_GENOMES); do f=$${F[$$g]}; s=$${S[$$g]}; \
		if [ ! -f "$$f" ]; then echo "  skip $$g ($$f not found)"; continue; fi; \
		for m in $(SWEEP_METHODS); do \
			echo ">> sweep $$m x $$g ($$s)  k: $(KSWEEP)"; out=$(OUT)/sweep/$$g.$$m.json; echo "[" > $$out; \
			for k in $(KSWEEP); do \
				case $$m in \
				minimizer) extra="$(SWEEP_W)";; \
				fracmh)    extra="$(SWEEP_MOD)";; \
				syncmer)   extra="$(SWEEP_S)";; \
				blend)     extra="$(SWEEP_W) $(SWEEP_B) $(SWEEP_N)";; \
				*)         extra="";; \
				esac; \
				$(BIN)/sketch-seeds.$$m "$$f" "$$s" $$k $$extra 2>/dev/null || \
				$(BIN)/sketch-seeds.$$m "$$f" $$k $$extra; \
			done >> $$out; \
			perl -0777 -i -pe 's/,\s*\z/\n/' $$out; \
			echo "]" >> $$out; echo "  -> $$out"; \
		done; \
	done; echo "-> $(OUT)/sweep/"

# plots
plot-uniqueness:
	@echo "Plotting k-mer uniqueness / gap ..."
	@ls $(OUT)/match-length/*.json >/dev/null 2>&1 || { echo "no match-length JSON (run: make matchlen)"; exit 1; }; \
	$(PYTHON) fit-uniqueness.py $(OUT)/match-length/*.json --out $(OUT)/kmer_uniqueness;
	@echo ""

sketch-plots:
	@echo "Plotting sketching uniqueness / gap ..."
	@declare -A S=( $(foreach g,$(GENOMES),[$(g)]="$($(g)_NAME)") ); \
	for g in $(GENOMES); do s=$${S[$$g]}; \
		echo ">> sketch-plot $$g ($$s)"; \
		$(PYTHON) sketch-plot.py "$$g" "$$s"; \
	done
	@echo ""

sketch-sweep-plots:
	@echo "Plotting sketch sweeps ..."
	@ls $(OUT)/sweep/*.json >/dev/null 2>&1 || { echo "no sweep JSON (run: make sweep)"; exit 1; }; \
	declare -A S=( $(foreach g,$(GENOMES),[$(g)]="$($(g)_NAME)") ); \
	for g in $(GENOMES); do s=$${S[$$g]}; \
		echo ">> sketch-sweep-plot $$g ($$s)"; \
		$(PYTHON) sketch-sweep-plot.py "$$g" "$$s"; \
	done
	@echo ""

plots: plot-uniqueness sketch-plots sketch-sweep-plots
	