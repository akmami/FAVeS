FASTA := ../data/chm13v2.0.fa

# Arguments
K_VALUES := 15 17 19 21 23 25 27 29 31
W_VALUES := 11 16 21 26 31 36 41 46 50

# Programs
# fa
BLEND := fuzzy-seeds
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
	@if [ -z "$(FASTA)" ]; then \
        echo "FASTA is not set or is empty. Exiting"; \
		exit 1; \
    fi
	@if [ ! -f $(FASTA) ]; then \
		echo "$(FASTA) does not exist. Exiting."; \
		exit 1; \
	fi

blend: 
	@mkdir -p ../lib
	gcc -O3 -mavx2 -msse4.1 -c ../sketch/blend.c -I blend
	ar rcs ../lib/libblend.a blend.o
	@rm -f blend.o

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
	@echo "$(FASTA) -> ../out/fa-blend/$(BLEND)-output.txt"
	@date > ../out/fa-blend/$(BLEND)-output.txt
	@echo "" >> ../out/fa-blend/$(BLEND)-output.txt

	@$(foreach k,$(K_VALUES), \
		$(foreach w,$(W_VALUES), \
			echo "Running k=$(k), w=$(w)" >> ../out/fa-blend/$(BLEND)-output.txt; \
			$(TIME) ../bin/$(BLEND) $(FASTA) $(k) $(w) >> ../out/fa-blend/$(BLEND)-output.txt 2>&1; \
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
	@echo "$(FASTA) -> ../out/fa-blend/$(BLEND)-output.txt"
	@date > ../out/fa-blend/$(BLEND)-output.txt
	@echo "" >> ../out/fa-blend/$(BLEND)-output.txt

	@echo "Running k=21, w=11, b=32" >> ../out/fa-blend/$(BLEND)-output.txt;
	$(TIME) ../bin/$(BLEND) $(FASTA) 21 11 32 >> ../out/fa-blend/$(BLEND)-output.txt 2>&1;
	@echo "" >> ../out/fa-blend/$(BLEND)-output.txt

	@echo "Running k=21, w=11, b=38" >> ../out/fa-blend/$(BLEND)-output.txt;
	$(TIME) ../bin/$(BLEND) $(FASTA) 21 11 38 >> ../out/fa-blend/$(BLEND)-output.txt 2>&1;
	@echo "" >> ../out/fa-blend/$(BLEND)-output.txt

	@echo "Running k=21, w=11, b=42" >> ../out/fa-blend/$(BLEND)-output.txt;
	$(TIME) ../bin/$(BLEND) $(FASTA) 21 11 42 >> ../out/fa-blend/$(BLEND)-output.txt 2>&1;
	@echo "" >> ../out/fa-blend/$(BLEND)-output.txt

	@echo "Running k=15, w=10, b=32" >> ../out/fa-blend/$(BLEND)-output.txt;
	$(TIME) ../bin/$(BLEND) $(FASTA) 15 10 32 >> ../out/fa-blend/$(BLEND)-output.txt 2>&1;
	@echo "" >> ../out/fa-blend/$(BLEND)-output.txt

	@echo "Running k=15, w=10, b=38" >> ../out/fa-blend/$(BLEND)-output.txt;
	$(TIME) ../bin/$(BLEND) $(FASTA) 15 10 38 >> ../out/fa-blend/$(BLEND)-output.txt 2>&1;
	@echo "" >> ../out/fa-blend/$(BLEND)-output.txt

	@echo "Running k=15, w=10, b=42" >> ../out/fa-blend/$(BLEND)-output.txt;
	$(TIME) ../bin/$(BLEND) $(FASTA) 15 10 42 >> ../out/fa-blend/$(BLEND)-output.txt 2>&1;
	@echo "" >> ../out/fa-blend/$(BLEND)-output.txt

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
	@rm -rf ../data/sim/*
	@echo "Output -> ../out/err-sim/$(SIMULATE)-output.txt"
	@date > $(OUT_DIR)/err-sim/$(SIMULATE)-output.txt
	@echo "" >> ../out/err-sim/$(SIMULATE)-output.txt

	$(TIME) ../bin/$(SIMULATE) -g -d ../data/sim >> ../out/err-sim/$(SIMULATE)-output.txt 2>&1
	@echo "" >> ../out/err-sim/$(SIMULATE)-output.txt

	$(TIME) ../bin/$(SIMULATE) -a --blend --short -d ../data/sim >> ../out/err-sim/$(SIMULATE)-output.txt 2>&1
	@echo "" >> ../out/err-sim/$(SIMULATE)-output.txt
	$(TIME) ../bin/$(SIMULATE) -a --blend -k 21 -w 11 -b 38 -d ../data/sim >> ../out/err-sim/$(SIMULATE)-output.txt 2>&1
	@echo "" >> ../out/err-sim/$(SIMULATE)-output.txt
	
	$(TIME) ../bin/$(SIMULATE) -a --blend --hifi -b 32 -d ../data/sim >> ../out/err-sim/$(SIMULATE)-output.txt 2>&1
	@echo "" >> ../out/err-sim/$(SIMULATE)-output.txt
	$(TIME) ../bin/$(SIMULATE) -a --blend --hifi -d ../data/sim >> ../out/err-sim/$(SIMULATE)-output.txt 2>&1
	@echo "" >> ../out/err-sim/$(SIMULATE)-output.txt

	$(TIME) ../bin/$(SIMULATE) -a --blend -k 15 -w 10 -b 32 -d ../data/sim >> ../out/err-sim/$(SIMULATE)-output.txt 2>&1
	@echo "" >> ../out/err-sim/$(SIMULATE)-output.txt
	$(TIME) ../bin/$(SIMULATE) -a --blend -k 15 -w 10 -b 38 -d ../data/sim >> ../out/err-sim/$(SIMULATE)-output.txt 2>&1
	@echo "" >> ../out/err-sim/$(SIMULATE)-output.txt
	
