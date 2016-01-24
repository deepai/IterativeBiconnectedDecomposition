IDIR =.
ODIR=obj
SRCDIR=.

EX_INCLUDE = ../include

CC=g++
CFLAGS=-I$(IDIR) -I$(EX_INCLUDE) -fopenmp -lgomp -g -w -std=c++11 -fpermissive

DFLAGS = -DMULTI_THREAD -DPRINT

#-DMULTI_THREAD for multi-thread support
#-DPRINT to print the bccs in a file

LIBS=-lm

_DEPS = *.h 
DEPS = $(patsubst %,$(EX_INCLUDE)/%,$(_DEPS))

_OBJ = helper.o mmio.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_SRC = *.cpp
SRC= $(patsubst %,$(SRCDIR)/%,$(_SRC))

_EX_SRC = mmio.cpp
EX_SRC = $(patsubst %,$(EX_INCLUDE)/%,$(_EX_SRC))

all: helper

$(ODIR)/%.o: $(SRC)
	$(CC) -c -o $@ $< $(CFLAGS) $(DFLAGS)

helper: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(DFLAGS)

$(ODIR)/mmio.o: $(EX_SRC) $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(DFLAGS)


.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
