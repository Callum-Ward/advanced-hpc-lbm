#Makefile

EXE=d2q9-bgk

CC=mpiicc
#CFLAGS= -std=c11 -Wall -O3
CFLAGS= -std=c11 -Wall -O3 -msse4 -mtune=native -march=native -funroll-loops --param max-unroll-times=4 -ffast-math
#CFLAGS= -O3 -msse4 -mtune=native -march=native -funroll-loops --param max-unroll-times=4 -ffast-math
#IFLAGS= -O3
IFLAGS= -O3 -xBROADWELL -restrict
LIBS = -lm

FINAL_STATE_FILE=./final_state.dat
AV_VELS_FILE=./av_vels.dat
REF_FINAL_STATE_FILE=check/128x128.final_state.dat
REF_AV_VELS_FILE=check/128x128.av_vels.dat

all: $(EXE)

$(EXE): $(EXE).c
	$(CC) $(IFLAGS) $^ $(LIBS) -o $@

test:
	mpiicc $(IFLAGS) test.c $(LIBS) -o test


check:
	python check/check.py --ref-av-vels-file=$(REF_AV_VELS_FILE) --ref-final-state-file=$(REF_FINAL_STATE_FILE) --av-vels-file=$(AV_VELS_FILE) --final-state-file=$(FINAL_STATE_FILE)
.PHONY: all check clean

clean:
	rm -f $(EXE)

cleant:
	rm -f test
