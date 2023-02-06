CC = c++
CFLAGS = -std=c++17 -Wpedantic -Wall -Wextra -Wshadow -Wnon-virtual-dtor \
	  -Wold-style-cast -Wcast-align -Wuseless-cast -Wsign-conversion \
	  -Wdouble-promotion -Wnull-dereference -Wmisleading-indentation \
	  -Wduplicated-cond -Wformat=2
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

PHONY += testLF
testLF:
	$(TARGET) LF

PHONY += testFORCE
testFORCE:
	$(TARGET) FORCE

PHONY += plot
FINDDATS = $(wildcard $(DAT)/*.dat)
RMSUFFIX := $(patsubst %.dat, %, $(FINDDATS))
RMDIR := $(patsubst ./dat/%, %, $(RMSUFFIX))
plot:
	for filename in $(RMDIR) ; do \
		gnuplot -e "filename='$$filename'" ./produceGif.plt ; \
	done

clean:
	rm $(OBJ)/*.o
	rm $(BIN)/*
	rm $(DAT)/*.dat
	rm *.gif

.PHONY: $(PHONY)
