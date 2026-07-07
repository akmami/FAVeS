CHM13_FASTA := ../data/chm13v2.0.fa
HG38_FASTA := ../data/human_v38.fasta

# Programs
# fa
SEEDS := sketch-seeds
SIMULATE := simulation

# Prefixes
BLEND := blend-seeds
MINIMIZER := minimizer-seeds
SYNCMER := syncmer-seeds
STROBEMER := strobemer-seeds
STROBEMER2 := strobemer2-seeds

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

syncmer: 
	@mkdir -p ../lib
	gcc -O3 -mavx2 -msse4.1 -c ../sketch/syncmer.c -I syncmer
	ar rcs ../lib/libsyncmer.a syncmer.o
	@rm -f syncmer.o

strobemer:
	@mkdir -p ../lib
	g++ -O3 -mavx2 -msse4.1 -std=c++17 -c ../sketch/strobemer/strobemer.cpp -I strobemer
	ar rcs ../lib/libstrobemer.a strobemer.o
	@rm -f strobemer.o

strobemer2:
	@mkdir -p ../lib
	g++ -O3 -mavx2 -msse4.1 -std=c++17 -c ../sketch/strobemer/strobemer2.cpp -I strobemer2
	ar rcs ../lib/libstrobemer2.a strobemer2.o
	@rm -f strobemer2.o

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
# BLEND
######################################################################
faBlend: check_fasta mkdir_bin mkdir_out blend
	$(GXX) $(CXXFLAGS) -mavx2 -msse4.1 -DBLEND -c $(SEEDS)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(SEEDS) $(SEEDS).o -L ../lib -lblend
	@rm $(SEEDS).o
	@mv $(SEEDS) $(EXECUTABLE_DIR)
	@mkdir -p ../out/fa-blend
	@echo "Output -> ../out/fa-blend/$(BLEND)-output.txt"
	@rm -f ../out/fa-blend/$(BLEND)-output.txt ../out/fa-blend/tmp.$(BLEND).*.log
	@echo "[" > ../out/fa-blend/$(BLEND)-output.txt
	@parallel -j $$(nproc) \
		'../bin/$(SEEDS) $(HG38_FASTA) {1} {2} {3} {4} > ../out/fa-blend/tmp.$(BLEND).{1}.{2}.{3}.{4}.log 2>&1' :::+ 21 15 :::+ 11 10 ::: 32 38 50 ::: 3 5 7
	@cat ../out/fa-blend/tmp.$(BLEND).*.log >> ../out/fa-blend/$(BLEND)-output.txt
	@echo "]" >> ../out/fa-blend/$(BLEND)-output.txt
	@rm -f ../out/fa-blend/tmp.$(BLEND).*.log

######################################################################
# MINIMIZERS
######################################################################
faMin: check_fasta mkdir_bin mkdir_out minimizer
	$(GXX) $(CXXFLAGS) -mavx2 -msse4.1 -DMINIMIZER -c $(SEEDS)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(SEEDS) $(SEEDS).o -L ../lib -lminimizer
	@rm $(SEEDS).o
	@mv $(SEEDS) $(EXECUTABLE_DIR)
	@mkdir -p ../out/fa-minimizer
	@echo "Output -> ../out/fa-minimizer/$(MINIMIZER)-output.txt"
	@rm -f ../out/fa-minimizer/$(MINIMIZER)-output.txt ../out/fa-minimizer/tmp.$(MINIMIZER).*.log
	@echo "[" > ../out/fa-minimizer/$(MINIMIZER)-output.txt
	@parallel -j 2 --colsep ' ' \
		'../bin/$(SEEDS) $(HG38_FASTA) {1} {2} > ../out/fa-minimizer/tmp.$(MINIMIZER).{1}.{2}.log 2>&1' ::: "21 11" "15 10"
	@cat ../out/fa-minimizer/tmp.$(MINIMIZER).*.log >> ../out/fa-minimizer/$(MINIMIZER)-output.txt
	@echo "]" >> ../out/fa-minimizer/$(MINIMIZER)-output.txt
	@rm -f ../out/fa-minimizer/tmp.$(MINIMIZER).*.log

######################################################################
# SYNCMERS
######################################################################
faSync: check_fasta mkdir_bin mkdir_out syncmer
	$(GXX) $(CXXFLAGS) -mavx2 -msse4.1 -DSYNCMER -c $(SEEDS)$(CXX) 
	$(GXX) $(CXXFLAGS) -o $(SEEDS) $(SEEDS).o -L ../lib -lsyncmer
	@rm $(SEEDS).o
	@mv $(SEEDS) $(EXECUTABLE_DIR)
	@mkdir -p ../out/fa-syncmer
	@echo "Output -> ../out/fa-syncmer/$(SYNCMER)-output.txt"
	@rm -f ../out/fa-syncmer/$(SYNCMER)-output.txt ../out/fa-syncmer/tmp.$(SYNCMER).*.log
	@echo "[" > ../out/fa-syncmer/$(SYNCMER)-output.txt
	@parallel -j 3 --colsep ' ' \
		'../bin/$(SEEDS) $(HG38_FASTA) {1} {2} > ../out/fa-syncmer/tmp.$(SYNCMER).{1}.{2}.log 2>&1' ::: "21 10" "15 10" "15 5"
	@cat ../out/fa-syncmer/tmp.$(SYNCMER).*.log >> ../out/fa-syncmer/$(SYNCMER)-output.txt
	@echo "]" >> ../out/fa-syncmer/$(SYNCMER)-output.txt
	@rm -f ../out/fa-syncmer/tmp.$(SYNCMER).*.log

######################################################################
# STROBEMERS
######################################################################
faStrobe: check_fasta mkdir_bin mkdir_out strobemer
	$(GXX) $(CXXFLAGS) -mavx2 -msse4.1 -DSTROBEMER -c $(SEEDS)$(CXX) 
	$(GXX) $(CXXFLAGS) -o $(SEEDS) $(SEEDS).o -L ../lib -lstrobemer
	@rm $(SEEDS).o
	@mv $(SEEDS) $(EXECUTABLE_DIR)
	@mkdir -p ../out/fa-strobemer
	@echo "Output -> ../out/fa-strobemer/$(STROBEMER)-output.txt"
	@rm -f ../out/fa-strobemer/$(STROBEMER)-output.txt ../out/fa-strobemer/tmp.$(STROBEMER).*.log
	@echo "[" > ../out/fa-strobemer/$(STROBEMER)-output.txt
	../bin/$(SEEDS) $(HG38_FASTA) 2 22 6 16 > ../out/fa-strobemer/$(STROBEMER)-output.txt 2>&1
	@echo "]" >> ../out/fa-strobemer/$(STROBEMER)-output.txt

######################################################################
# STROBEMERS2
######################################################################
faStrobe2: check_fasta mkdir_bin mkdir_out strobemer2
	$(GXX) $(CXXFLAGS) -mavx2 -msse4.1 -DSTROBEMER2 -c $(SEEDS)$(CXX) 
	$(GXX) $(CXXFLAGS) -o $(SEEDS) $(SEEDS).o -L ../lib -lstrobemer2
	@rm $(SEEDS).o
	@mv $(SEEDS) $(EXECUTABLE_DIR)
	@mkdir -p ../out/fa-strobemer2
	@echo "Output -> ../out/fa-strobemer2/$(STROBEMER2)-output.txt"
	@rm -f ../out/fa-strobemer2/$(STROBEMER2)-output.txt ../out/fa-strobemer2/tmp.$(STROBEMER2).*.log
	@echo "[" > ../out/fa-strobemer2/$(STROBEMER2)-output.txt
	../bin/$(SEEDS) $(HG38_FASTA) 250 17 > ../out/fa-strobemer2/$(STROBEMER2)-output.txt 2>&1
	@echo "]" >> ../out/fa-strobemer2/$(STROBEMER2)-output.txt


######################################################################
# BLEND ERROR TOLERANCE EXPERIMENT
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
	@parallel -j $$(nproc) \
		'../bin/$(SIMULATE) {1} {2} {3} {4} > ../out/err-sim/tmp.$(SIMULATE).{1}.{2}.{3}.{4}.log 2>&1' :::+ 21 15 :::+ 11 10 ::: 32 38 50 ::: 3 5 7
	@cat ../out/err-sim/tmp.$(SIMULATE).*.log >> ../out/err-sim/$(SIMULATE)-output.txt
	@echo "]" >> ../out/err-sim/$(SIMULATE)-output.txt
	@rm -f ../out/err-sim/tmp.$(SIMULATE).*.log
