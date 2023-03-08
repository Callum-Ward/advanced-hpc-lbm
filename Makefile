#Makefile

EXE=d2q9-bgk

CC=icc

IFLAGS= -O3 -xBROADWELL -Ofast -restrict
LIBS = -lm
ILIBS = -lm -qopenmp

FINAL_STATE_FILE=./final_state.dat
AV_VELS_FILE=./av_vels.dat
REF_FINAL_STATE_FILE=check/128x256.final_state.dat
REF_AV_VELS_FILE=check/128x256.av_vels.dat

all: $(EXE)

$(EXE): $(EXE).c
	$(CC) $(IFLAGS) $^ $(ILIBS) -o $@



check:
	python check/check.py --ref-av-vels-file=$(REF_AV_VELS_FILE) --ref-final-state-file=$(REF_FINAL_STATE_FILE) --av-vels-file=$(AV_VELS_FILE) --final-state-file=$(FINAL_STATE_FILE)
.PHONY: all check clean

clean:
	rm -f $(EXE)
