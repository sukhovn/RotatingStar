.PHONY: all clean

DOMAINS_DIR = Domains
ROUTINES_DIR = Routines
COLLOCATION_DIR = Collocation
STAR_CODE_DIR = Star
# EXTRAPOLATION_DIR = Extrapolation

DOMAINS_SRC := $(wildcard $(DOMAINS_DIR)/*.cpp)
DOMAINS_BINS := $(patsubst $(DOMAINS_DIR)/%.cpp, obj/domains/%.o, $(DOMAINS_SRC))

ROUTINES_SRC := $(wildcard $(ROUTINES_DIR)/*.cpp)
ROUTINES_BINS := $(patsubst $(ROUTINES_DIR)/%.cpp, obj/routines/%.o, $(ROUTINES_SRC))

COLLOCATION_SRC := $(wildcard $(COLLOCATION_DIR)/*.cpp)
COLLOCATION_BINS := $(patsubst $(COLLOCATION_DIR)/%.cpp, obj/colloc/%.o, $(COLLOCATION_SRC))

STAR_CODE_SRC := $(wildcard $(STAR_CODE_DIR)/*.cpp)
STAR_CODE_BINS := $(patsubst $(STAR_CODE_DIR)/%.cpp, obj/star/%.o, $(STAR_CODE_SRC))

# EXTRAPOLATION_SRC := $(wildcard $(EXTRAPOLATION_DIR)/*.cpp)
# EXTRAPOLATION_BINS := $(patsubst $(EXTRAPOLATION_DIR)/%.cpp, obj/extra/%.o, $(EXTRAPOLATION_SRC))

HEAD_FOLDERS = $(DOMAINS_DIR) $(ROUTINES_DIR) $(COLLOCATION_DIR) $(STAR_CODE_DIR) # $(EXTRAPOLATION_DIR)
HEAD_LINK = $(patsubst %, -I%, $(HEAD_FOLDERS))

CC = g++ # compiler to use
CFLAGS = -c -O3 -Iobj $(HEAD_LINK) # flags to use at the compliation
LINKERFLAG = -lm -lfftw3 -fopenmp -lblas -llapack -llapacke
BINARIES = $(DOMAINS_BINS) $(ROUTINES_BINS) $(COLLOCATION_BINS) $(STAR_CODE_BINS) # $(EXTRAPOLATION_BINS)

all: prog_initial prog_multidomain # prog_extrapolator

obj:
	@echo "Making directory for binaries..."
	mkdir obj
	mkdir obj/domains
	mkdir obj/routines
	mkdir obj/colloc
	mkdir obj/star
	mkdir obj/extra
	xxd -i ./$(STAR_CODE_DIR)/solution_parameters.dat ./obj/solution_parameters.h
	xxd -i ./$(STAR_CODE_DIR)/init_setup_parameters.dat ./obj/init_setup_parameters.h
# 	xxd -i ./$(EXTRAPOLATION_DIR)/extrapolation_save_par.dat ./obj/extrapolation_save_par.h

obj/domains/%.o: $(DOMAINS_DIR)/%.cpp obj
	@echo "Compiling domains binary" $@ "..."
	${CC} ${CFLAGS} -o $@ $<

obj/routines/%.o: $(ROUTINES_DIR)/%.cpp obj
	@echo "Compiling routine binary" $@ "..."
	${CC} ${CFLAGS} -o $@ $<

obj/colloc/%.o: $(COLLOCATION_DIR)/%.cpp obj
	@echo "Compiling collocation binary" $@ "..."
	${CC} ${CFLAGS} -o $@ $<

obj/star/%.o: $(STAR_CODE_DIR)/%.cpp obj
	@echo "Compiling star binary" $@ "..."
	${CC} ${CFLAGS} -o $@ $<

# obj/extra/%.o: $(EXTRAPOLATION_DIR)/%.cpp obj
# 	@echo "Compiling extrapolation binary" $@ "..."
# 	${CC} ${CFLAGS} -o $@ $<

obj/initial.o: initial.cpp obj
	@echo "Compiling the initial guess solver..."
	${CC} ${CFLAGS} -o obj/initial.o initial.cpp

obj/multidomain.o: multidomain.cpp obj
	@echo "Compiling the multidomain solver..."
	${CC} ${CFLAGS} -o obj/multidomain.o multidomain.cpp

# obj/extrapolator.o: extrapolator.cpp obj
# 	@echo "Compiling the extrapolator..."
# 	${CC} ${CFLAGS} -o obj/extrapolator.o extrapolator.cpp

prog_initial: $(BINARIES) obj/initial.o
	@echo "Linking the initial guess solver..."
	${CC} -o initial $(BINARIES) obj/initial.o $(LINKERFLAG)

prog_multidomain: $(BINARIES) obj/multidomain.o
	@echo "Linking the multidomain solver..."
	${CC} -o multidomain $(BINARIES) obj/multidomain.o $(LINKERFLAG)

# prog_extrapolator: $(BINARIES) obj/extrapolator.o
# 	@echo "Linking the extrapolator..."
# 	${CC} -o extrapolator $(BINARIES) obj/extrapolator.o $(LINKERFLAG)

clean:
	@echo "Cleaning up..."
	rm -rvf obj
	rm -vf initial
	rm -vf multidomain
	rm -vf extrapolator