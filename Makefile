CC = c++
CFLAGS = -std=c++17 -Wpedantic -Wall -Wextra -Wnon-virtual-dtor \
	  -Wold-style-cast -Wcast-align -Wuseless-cast \
	  -Wdouble-promotion -Wmisleading-indentation \
	  -Wduplicated-cond -Wformat=2 -O2
LIBFLAGS = -lconfig++
DBGFLAGS = -g -O0 -fsanitize=address -fsanitize=bounds -lubsan
OPTFLAG = -O2

BIN = ./bin
SRC = ./src
OBJ = ./obj
SET = ./stg
DAT = ./dat

TARGET = $(BIN)/fvSim

INCS = $(wildcard $(SRC)/*.H)
SRCS = $(wildcard $(SRC)/*.C)
OBJS = $(patsubst $(SRC)/%.C, $(OBJ)/%.o, $(SRCS))
INCDIRS = -I./ $(addprefix -I, $(SRC))

PHONY := $(TARGET)
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBFLAGS)

$(OBJ)/%.o: $(SRC)/%.C $(INCS)
	$(CC) $(CFLAGS) -o $@ -c $< $(INC_DIRS)

PHONY += debug
debug: 
	CFLAGS += $(DBGFLAGS)
	$(TARGET)

PHONY += release
release: 
	CFLAGS += $(OPTFLAG)
	$(TARGET)

PHONY += refactor
refactor: 
	CFLAGS += $(REFFLAGS)
	$(TARGET)

PHONY += runall
runall:
	$(TARGET) HLLC 1

PHONY += plotall
plotall:
	gnuplot ./plt/plotAllMHDTests.plt

PHONY += LF
LF:
	$(TARGET) LF

PHONY += FORCE
FORCE:
	$(TARGET) FORCE

PHONY += SLIC
SLIC:
	$(TARGET) SLIC

PHONY += HLL
HLL:
	$(TARGET) HLL

PHONY += HLLC
HLLC:
	$(TARGET) HLLC

PHONY += plotGifs
FINDDATS = $(wildcard $(DAT)/*.dat)
RMSUFFIX := $(patsubst %.dat, %, $(FINDDATS))
RMDIR := $(patsubst ./dat/%, %, $(RMSUFFIX))
plotGifs:
	for filename in $(RMDIR) ; do \
		gnuplot -e "filename='$$filename'" ./plt/produceGif.plt ; \
	done

PHONY += plotGifs1D
FINDDATS = $(wildcard $(DAT)/*1D.dat)
RMSUFFIX := $(patsubst %.dat, %, $(FINDDATS))
RMDIR := $(patsubst ./dat/%, %, $(RMSUFFIX))
plotGifs1D:
	for filename in $(RMDIR) ; do \
		gnuplot -e "filename='$$filename'" ./plt/produceGif1D.plt ; \
	done

PHONY += plotF
FINDDATS = $(wildcard $(DAT)/*.dat)
RMSUFFIX := $(patsubst %.dat, %, $(FINDDATS))
RMDIR := $(patsubst ./dat/%, %, $(RMSUFFIX))
plotF:
	for filename in $(RMDIR) ; do \
		gnuplot -e "filename='$$filename'" ./plt/plotFinal.plt ; \
	done

PHONY += plot1D
FINDDATS = $(wildcard $(DAT)/*1D.dat)
RMSUFFIX := $(patsubst %.dat, %, $(FINDDATS))
RMDIR := $(patsubst ./dat/%, %, $(RMSUFFIX))
plot1D:
	for filename in $(RMDIR) ; do \
		gnuplot -e "filename='$$filename'" ./plt/plotFinal1D.plt ; \
	done

PHONY += clean
clean:
	rm $(OBJ)/*.o
	rm $(BIN)/*

PHONY += cleanresults
cleanresults:
	rm $(DAT)/*.dat
	rm plots/*.gif
	rm plots/*.png

.PHONY: $(PHONY)
