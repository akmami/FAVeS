TARGET := faves
SRCS := $(wildcard *.c)
OBJS := $(SRCS:.c=.o)
CURRENT_DIR := $(shell pwd)

CC := gcc
WFA2_DIR := WFA2-lib
CFLAGS := -O3 -mavx2 -msse4.1 -Wall -Wextra -Wpedantic -I $(WFA2_DIR)
LDFLAGS := -lm -pthread -lz

# ========================================
#  Build Rules
# ========================================

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -L lib -lblend -L $(WFA2_DIR)/lib -lwfa $(LDFLAGS)
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
	$(CC) $(CFLAGS) -c sketch/blend.c -I blend
	ar rcs lib/libblend.a blend.o
	@rm -f blend.o

wfa:
	cd $(WFA2_DIR) && make && cd ..

install: clean wfa blend $(TARGET)
