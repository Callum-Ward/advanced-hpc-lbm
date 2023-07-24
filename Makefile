#Makefile
# std::execution port: https://github.com/nvidia/stdexec
# std::mdspan port: https://github.com/kokkos/mdspan

EXE=d2q9-bgk
EXETEST=d2q9-bgk-test

CC=icc
IFLAGS= -O3 -xBROADWELL -Ofast -restrict
LIBS = -lm
ILIBS = -lm -qopenmp

GC=g++
CFLAGS= -std=c++20 -Ofast -march=native -DNDEBUG
CLIBS= -ltbb

NC=nvc++
NFLAGS= -std=c++20 -O4 -fast -march=native -DNDEBUG -Mllvm-fast

FINAL_STATE_FILE=./final_state.dat
AV_VELS_FILE=./av_vels.dat
REF_FINAL_STATE_FILE=check/128x256.final_state.dat
REF_AV_VELS_FILE=check/128x256.av_vels.dat

all: $(EXE)

$(EXE): $(EXE).cpp
	$(GC) $(CFLAGS) $^ $(CLIBS) -o $@

test: $(EXETEST)

$(EXETEST): $(EXETEST).cpp
	$(GC) $(CFLAGS) $^ $(CLIBS) -o $@

d2q9-bgk-nvc: $(EXE).cpp
	$(NC) $(NFLAGS) $^ -o $@

d2q9-bgk-icc: $(EXE).c
	$(CC) $(IFLAGS) $^ $(ILIBS) -o $@




check:
	python check/check.py --ref-av-vels-file=$(REF_AV_VELS_FILE) --ref-final-state-file=$(REF_FINAL_STATE_FILE) --av-vels-file=$(AV_VELS_FILE) --final-state-file=$(FINAL_STATE_FILE)

clean:
	rm -f $(EXE) d2q9-bgk-nvc d2q9-bgk-icc	

.PHONY: all check clean


