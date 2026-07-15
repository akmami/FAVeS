TARGET := faves
SRCS := $(wildcard *.c)
OBJS := $(SRCS:.c=.o)
CURRENT_DIR := $(shell pwd)

CC := gcc
CXX := g++
WFA2_DIR := WFA2-lib
CFLAGS := -O3 -mavx2 -msse4.1 -Wall -Wextra -Wpedantic -I $(WFA2_DIR)
CXXFLAGS := -O3 -mavx2 -msse4.1 -std=c++17
LDFLAGS := -lm -pthread -lz

SKETCH_LIBS := -L lib -lblend -lminimizer -lsyncmer -llcp
CXX_RUNTIME := -lstdc++

# ========================================
#  Build Rules
# ========================================

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(SKETCH_LIBS) -L $(WFA2_DIR)/lib -lwfa $(LDFLAGS) $(CXX_RUNTIME)
	rm -f $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# ========================================
#  Utility Targets
# ========================================

clean:
	@echo "Cleaning"
	rm -f $(OBJS)
	rm -f $(TARGET)
	cd $(WFA2_DIR) && make clean && cd ..

blend:
	@mkdir -p lib
	$(CC) $(CFLAGS) -c sketch/blend.c -o blend.o
	ar rcs lib/libblend.a blend.o
	@rm -f blend.o

minimizer:
	@mkdir -p lib
	$(CC) $(CFLAGS) -c sketch/minimizer.c -o minimizer.o
	ar rcs lib/libminimizer.a minimizer.o
	@rm -f minimizer.o

syncmer:
	@mkdir -p lib
	$(CC) $(CFLAGS) -c sketch/syncmer.c -o syncmer.o
	ar rcs lib/libsyncmer.a syncmer.o
	@rm -f syncmer.o

lcp:
	@mkdir -p lib
	$(CXX) $(CXXFLAGS) -c sketch/lps.c -o lps.o
	ar rcs lib/liblcp.a lps.o
	@rm -f lps.o

sketches: blend minimizer syncmer lcp

wfa:
	cd $(WFA2_DIR) && make && cd ..

install: clean wfa sketches $(TARGET)