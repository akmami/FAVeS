CHM13_FASTA := ../data/chm13v2.0.fa
HG38_FASTA := ../data/human_v38.fasta

# Arguments
K_VALUES := 15 17 19 21 23 25 27 29 31
W_VALUES := 11 16 21 26 31 36 41 46 50

# Programs
# fa
BLEND := fuzzy-seeds
MINIMIZER := min-seeds
SIMULATE := simulation

# Directories
CURRENT_DIR := $(shell pwd)
DATA_DIR := $(CURRENT_DIR)/../data
EXECUTABLE_DIR := $(CURRENT_DIR)/../bin
OUT_DIR := $(CURRENT_DIR)/../out

# Extension
CXX := .cpp

# Compiler
GXX := g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -Wpedantic
BLENDLIB = -L ../lib -lblend 
MINIMIZERLIB = -L ../lib -lminimizer
TIME := /usr/bin/time -v

.PHONY: all %

%:
	@:

all: mkdir_bin mkdir_out

mkdir_bin:
	@if [ ! -d "$(EXECUTABLE_DIR)" ]; then \
		echo "$(EXECUTABLE_DIR) does not exist. Creating..."; \
		mkdir $(EXECUTABLE_DIR); \
	fi

mkdir_out:
	@if [ ! -d "$(OUT_DIR)" ]; then \
		echo "$(OUT_DIR) does not exist. Creating..."; \
		mkdir $(OUT_DIR); \
	fi

check_fasta:
	@if [ -z "$(CHM13_FASTA)" ]; then \
        echo "FASTA is not set or is empty. Exiting"; \
		exit 1; \
    fi
	@if [ ! -f $(CHM13_FASTA) ]; then \
		echo "$(CHM13_FASTA) does not exist. Exiting."; \
		exit 1; \
	fi

blend: 
	@mkdir -p ../lib
	gcc -O3 -mavx2 -msse4.1 -c ../sketch/blend.c -I blend
	ar rcs ../lib/libblend.a blend.o
	@rm -f blend.o

minimizer: 
	@mkdir -p ../lib
	gcc -O3 -mavx2 -msse4.1 -c ../sketch/minimizer.c -I minimizer
	ar rcs ../lib/libminimizer.a minimizer.o
	@rm -f minimizer.o

clean:
	@if [ -d "$(OUT_DIR)" ]; then \
		echo "Removing $(OUT_DIR)"; \
		rm -r $(OUT_DIR); \
	fi
	@if [ -d "$(EXECUTABLE_DIR)" ]; then \
		echo "Removing $(EXECUTABLE_DIR)"; \
		rm -r $(EXECUTABLE_DIR); \
	fi
	@echo "Cleaned."

######################################################################
# FUZZY SEEDS (BLEND-ALL)
######################################################################
faBlendAll: check_fasta mkdir_bin mkdir_out blend
	$(GXX) $(CXXFLAGS) -mavx2 -msse4.1 -c $(BLEND)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(BLEND) $(BLEND).o $(BLENDLIB)
	@rm $(BLEND).o
	@mv $(BLEND) $(EXECUTABLE_DIR)

	@mkdir -p ../out/fa-blend
	@echo "$(CHM13_FASTA) -> ../out/fa-blend/$(BLEND)-output.txt"
	@date > ../out/fa-blend/$(BLEND)-output.txt
	@echo "" >> ../out/fa-blend/$(BLEND)-output.txt

	@$(foreach k,$(K_VALUES), \
		$(foreach w,$(W_VALUES), \
			echo "Running k=$(k), w=$(w)" >> ../out/fa-blend/$(BLEND)-output.txt; \
			$(TIME) ../bin/$(BLEND) $(CHM13_FASTA) $(k) $(w) >> ../out/fa-blend/$(BLEND)-output.txt 2>&1; \
			echo "" >> ../out/fa-blend/$(BLEND)-output.txt; \
		) \
	)

######################################################################
# FUZZY SEEDS (BLEND-SINGLE)
######################################################################
faBlend: check_fasta mkdir_bin mkdir_out blend
	$(GXX) $(CXXFLAGS) -mavx2 -msse4.1 -c $(BLEND)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(BLEND) $(BLEND).o $(BLENDLIB)
	@rm $(BLEND).o
	@mv $(BLEND) $(EXECUTABLE_DIR)
	@mkdir -p ../out/fa-blend
	@echo "Output -> ../out/fa-blend/$(BLEND)-output.txt"
	@rm -f ../out/fa-blend/$(BLEND)-output.txt ../out/fa-blend/tmp.$(BLEND).*.log
	@echo "[" > ../out/fa-blend/$(BLEND)-output.txt
	@printf "%s\n" \
		"21 11 32 3" "21 11 32 5" "21 11 32 7" \
		"21 11 37 3" "21 11 37 5" "21 11 37 7" \
		"21 11 42 3" "21 11 42 5" "21 11 42 7" \
		"15 10 32 3" "15 10 32 5" "15 10 32 7" \
		"15 10 37 3" "15 10 37 5" "15 10 37 7" \
		"15 10 42 3" "15 10 42 5" "15 10 42 7" \
	| xargs -n 4 -P $$(nproc) bash -c '\
		a=$$1; b=$$2; c=$$3; d=$$4; \
		f="../out/fa-blend/tmp.$(BLEND).$${a}.$${b}.$${c}.$${d}.log"; \
		{ \
		  ../bin/$(BLEND) $(HG38_FASTA) $$a $$b $$c $$d; \
		} >"$$f" 2>&1 \
	' _ 
	@cat ../out/fa-blend/tmp.$(BLEND).*.log >> ../out/fa-blend/$(BLEND)-output.txt
	@echo "]" >> ../out/fa-blend/$(BLEND)-output.txt
	@rm -f ../out/fa-blend/tmp.$(BLEND).*.log
# 	@../bin/$(BLEND) $(HG38_FASTA) 21 11 32 7 yes
# 	@mkdir -p ../out/fa-blend-span
# 	@mv ../out/relaxed_span_*.bed ../out/fa-blend-span

######################################################################
# ERROR TOLERANCE EXPERIMENT
######################################################################
errSim: mkdir_bin mkdir_out blend
	$(GXX) $(CXXFLAGS) -mavx2 -msse4.1 -c $(SIMULATE)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(SIMULATE) $(SIMULATE).o $(BLENDLIB)
	@rm $(SIMULATE).o
	@mv $(SIMULATE) $(EXECUTABLE_DIR)
	@mkdir -p ../out/err-sim
	@mkdir -p ../data/sim
	@echo "Output -> ../out/err-sim/$(SIMULATE)-output.txt"
	@rm -f ../out/err-sim/$(SIMULATE)-output.txt ../out/err-sim/tmp.$(SIMULATE).*.log
	@echo "[" > ../out/err-sim/$(SIMULATE)-output.txt
	@printf "%s\n" \
		"21 11 32 3" "21 11 32 5" "21 11 32 7" \
		"21 11 37 3" "21 11 37 5" "21 11 37 7" \
		"21 11 42 3" "21 11 42 5" "21 11 42 7" \
		"15 10 32 3" "15 10 32 5" "15 10 32 7" \
		"15 10 37 3" "15 10 37 5" "15 10 37 7" \
		"15 10 42 3" "15 10 42 5" "15 10 42 7" \
	| xargs -n 4 -P $$(nproc) bash -c '\
		a=$$1; b=$$2; c=$$3; d=$$4; \
		f="../out/err-sim/tmp.$(SIMULATE).$${a}.$${b}.$${c}.$${d}.log"; \
		{ \
		  ../bin/$(SIMULATE) $$a $$b $$c $$d; \
		} >"$$f" 2>&1 \
	' _ 
	@cat ../out/err-sim/tmp.$(SIMULATE).*.log >> ../out/err-sim/$(SIMULATE)-output.txt
	@echo "]" >> ../out/err-sim/$(SIMULATE)-output.txt
	@rm -f ../out/err-sim/tmp.$(SIMULATE).*.log


######################################################################
# MINIMIZERS
######################################################################
faMin: check_fasta mkdir_bin mkdir_out minimizer
	$(GXX) $(CXXFLAGS) -mavx2 -msse4.1 -c $(MINIMIZER)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(MINIMIZER) $(MINIMIZER).o $(MINIMIZERLIB)
	@rm $(MINIMIZER).o
	@mv $(MINIMIZER) $(EXECUTABLE_DIR)
	@mkdir -p ../out/fa-minimizer
	@mkdir -p ../out/fa-minimizer-span
	@echo "Output -> ../out/fa-minimizer/$(MINIMIZER)-output.txt"
	@rm -f ../out/fa-minimizer/$(MINIMIZER)-output.txt ../out/fa-minimizer/tmp.$(MINIMIZER).*.log
	@echo "[" > ../out/fa-minimizer/$(MINIMIZER)-output.txt
	@printf "%s\n" \
		"21 11" "15 10" \
	| xargs -n 2 -P $$(nproc) bash -c '\
		a=$$1; b=$$2; \
		f="../out/fa-minimizer/tmp.$(MINIMIZER).$${a}.$${b}.log"; \
		{ \
		  ../bin/$(MINIMIZER) $(HG38_FASTA) $$a $$b yes; \
		} >"$$f" 2>&1 \
	' _ 
	@cat ../out/fa-minimizer/tmp.$(MINIMIZER).*.log >> ../out/fa-minimizer/$(MINIMIZER)-output.txt
	@echo "]" >> ../out/fa-minimizer/$(MINIMIZER)-output.txt
	@rm -f ../out/fa-minimizer/tmp.$(MINIMIZER).*.log