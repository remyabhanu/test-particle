# Master Makefile

# Declare DaViTpy working directory
DAVITDIR = $(DAVITPY)

# Declare local compilers, or use system defaults
FC   = gfortran
F77  = gfortran
MPI  = mpif90
F2PY = f2py

# If you have problems compiling certain routines, your system may not support
# some of these flags.  However, it is not recommended to change them.
FC_FLAGS = -O2 -fPIC
F77_FLAGS = $(FC_FLAGS) -fbacktrace -fno-automatic
MPI_FLAGS = -O2 -fbacktrace -fno-automatic

# Declare directories which contain C or Fortran code that must be compiled
# NOTE: IRGR and MSIS are currently pre-compiled and so not included here
TSYDIR = tsyganenko

all: clean build

build:
	(cd $(TSYDIR); make F77=$(F77) F2PY=$(F2PY) OPT_FLAGS="$(F77_FLAGS)")

clean:
	(cd $(TSYDIR); make clean)
test:
	(cd $(TSYDIR); make test)
test_clean:
	(cd $(TSYDIR); make test_clean)
