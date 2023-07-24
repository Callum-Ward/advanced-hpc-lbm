#Makefile
# std::execution port: https://github.com/nvidia/stdexec
# std::mdspan port: https://github.com/kokkos/mdspan

EXE=d2q9-bgk

CC=icc
IFLAGS= -O3 -xBROADWELL -Ofast -restrict
LIBS = -lm
ILIBS = -lm -qopenmp

GC=g++
CFLAGS= -std=c++20 -Ofast -march=native -DNDEBUG
CLIBS= -ltbb

NC=nvc++
NFLAGS= -stdpar=multicore -std=c++20 -O4 -fast -march=native -Mllvm-fast -DNDEBUG

FINAL_STATE_FILE=./final_state.dat
AV_VELS_FILE=./av_vels.dat
REF_FINAL_STATE_FILE=check/128x128.final_state.dat
REF_AV_VELS_FILE=check/128x128.av_vels.dat

all: $(EXE)

$(EXE): $(EXE).cpp
	$(NC) $(NFLAGS) $^ -o $(EXE)

g++: $(EXE).cpp
	$(GC) $(CFLAGS) $^ $(CLIBS) -o $(EXE)

icc: $(EXE).c
	$(CC) $(IFLAGS) $^ $(ILIBS) -o $(EXE)

run:
	./$(EXE) input_128x128.params obstacles_128x128.dat


check:
	python check/check.py --ref-av-vels-file=$(REF_AV_VELS_FILE) --ref-final-state-file=$(REF_FINAL_STATE_FILE) --av-vels-file=$(AV_VELS_FILE) --final-state-file=$(FINAL_STATE_FILE)

clean:
	rm -f $(EXE)

.PHONY: all check clean


