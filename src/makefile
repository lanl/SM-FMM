# Compiler
FC = gfortran

# Source files
SRC = fftw3.f90 \
  fourier_functions_r2c.f90 \
  global.f90 \
  types.f90 \
  tensor_functions.f90 \
  various_functions.f90 \
  IO_functions.f90 \
  main.f90
OBJ = $(SRC:.f90=.o)

# LAPACK location. This line may be left blank on some systems.
LAPACK_Windows := C:\Users\304378\.conda\envs\LSEVPFFT\Library\bin
LAPACK_Darwin := 
LAPACK_Linux := 

# FFTW location
FFTW3_Windows := C:\Users\304378\.conda\envs\LSEVPFFT\Library\lib
FFTW3_Darwin := /Users/miroslavzecevic/miniconda3/envs/fortran/
FFTW3_Linux := /home/miroslav/fftw3

# Detect operating system
ifeq ($(OS),Windows_NT)
  detected_OS := Windows
else
  detected_OS := $(shell sh -c 'uname 2>/dev/null || echo Unknown')
endif

OS:
	@echo $(detected_OS)

# Set folders
ifeq ($(detected_OS),Windows) 
  LAPACK := $(LAPACK_Windows)
  FFTW3 := $(FFTW3_Windows)
endif
ifeq ($(detected_OS),Darwin) 
  LAPACK := $(LAPACK_Darwin)
  FFTW3 := $(FFTW3_Darwin)
endif
ifeq ($(detected_OS),Linux) 
  LAPACK := $(LAPACK_Linux)
  FFTW3 := $(FFTW3_Linux)
endif

# Check for fftw folder
ifeq ("$(wildcard $(FFTW3))","")
  $(error Directory "$(FFTW3)" does not exist, check FFTW3 flag in makefile)
endif

# Check for fftw3.f03
ifeq ("$(wildcard fftw3.f03)","")
  $(error Missing file "fftw3.f03" - copy the file to this directory from your fftw download (/fftw-3.x.x/api/fftw3.f03))
endif

# Executable name
ifeq ($(detected_OS),Windows)
  EXEC = FMM.exe
else
  EXEC = FMM.out
endif

# openmp make option
OMPFLAGS := -fopenmp 
ifeq ($(detected_OS),Windows)
	FFTWOMP :=
else
	FFTWOMP := -lfftw3_threads -lm
endif

# Setup flags
ifeq ($(strip $(LAPACK)),)
  LAPACKFLAGS := -llapack
else
  LAPACKFLAGS := -L$(LAPACK) -llapack
endif
ifeq ($(detected_OS),Windows)
  FFTWFLAGS := -L${FFTW3} -lfftw3 -lm $(FFTWOMP)
else
  FFTWFLAGS := -L${FFTW3}/lib -lfftw3 -lm $(FFTWOMP)
endif

# Help
help: 
	@echo Type make debug for a debug build with runtime checks.
	@echo Type make opt for a fast, optimized build.
	@echo Type make openmp for a fast, optimized build with openmp parallelization.
	@echo Type make clean to remove all build products.

# Flags
FFLAGS_BASIC := -m64 -fdefault-real-8 -fdefault-double-8 -fno-align-commons
FFLAGS_LIB := $(LAPACKFLAGS) $(FFTWFLAGS)
FFLAGS_DEBUG := -Wall -Wline-truncation -Wextra -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface \
            -Wunused-parameter -Wconversion -fwhole-file -pedantic -fbacktrace -fbounds-check -fcheck=all -ftrapv -O0 -C -g -fimplicit-none \
			-ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan -ffree-line-length-0
FFLAGS_OPT := -O3 -ffast-math -fno-protect-parens -fno-signed-zeros -fvariable-expansion-in-unroller -funroll-loops

# Default target
all: debug

# Debug build
debug: FFLAGS = $(FFLAGS_BASIC) $(FFLAGS_DEBUG) $(FFLAGS_LIB) 
debug: FFLAGS_LINK = $(FFLAGS_LIB) 
debug: FFLAGS_LINK_BEFORE := 
debug: $(EXEC)

# Optimized build
opt: FFLAGS = $(FFLAGS_BASIC) $(FFLAGS_OPT) $(FFLAGS_LIB) 
opt: FFLAGS_LINK = $(FFLAGS_LIB) 
opt: FFLAGS_LINK_BEFORE := 
opt: $(EXEC)

# Parallel build
openmp: FFLAGS = $(FFLAGS_BASIC) $(FFLAGS_OPT) $(FFLAGS_LIB) $(OMPFLAGS)
openmp: FFLAGS_LINK = $(FFLAGS_LIB)
openmp: FFLAGS_LINK_BEFORE = $(OMPFLAGS)
openmp: $(EXEC)

# Linking
$(EXEC): $(OBJ)
	$(FC) $(FFLAGS_LINK_BEFORE) -o $@ $^ $(FFLAGS_LINK)

# Compilation rule
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f *.o *.mod *.smod $(EXEC)

