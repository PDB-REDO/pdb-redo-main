firstTarget: all


BIN = tools
SRC = src

SOURCES = $(wildcard $(SRC)/*.f)
SOURCES := $(SOURCES:$(SRC)/%.f=$(SRC)/%.f)
EXE = $(SOURCES:$(SRC)/%.f=$(BIN)/%)


FC = gfortran
FFLAGS =


.PHONY: all clean

all: $(EXE)

$(BIN)/%: $(SRC)/%.f
	$(FC) $< $(FFLAGS) -o $@

clean:
	rm -f $(EXE)
