#------------------------------------------------------------------------------#
# Makefile                                                                     #
#                                                                              #
# This Makefile is structured fairly linearly -- at the top are some reasonable#
# default settings.  Following, we override these defaults and set up specific #
# configurations for the various machines we use.  At the bottom are the nuts  #
# and bolts of the compilation which likely needn't be touched frequently.     #
#                                                                              #
# TARGETS:                                                                     #
# nompi    => pnfam_nompi.x                                                    #
# mpi      => pnfam_mpi.x                                                      #
# beta     => betadecay.x                                                      #
# hfbinfo  => hfbinfo.x                                                        #
# tests    => tests/                                                           #
#                                                                              #
# M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15                        #
# Modified by T. Li, Michigan State Univ, 2020
#------------------------------------------------------------------------------#

#==============================================================================#
#                           Compile-Time Options                               #
#==============================================================================#
COMPILER = ifort 
USE_ADLB = 0
USE_MPI  = 1

# Use the hamiltonian module optimized with BLAS or not
USE_HBLAS = 1

# Can be overridden in the make command, e.g. "make openmp=1 all"
VERBOSE   = 0
debug     = 0
openmp    = 0

# Locations for dependencies:
# Basic Linear Algebra Subroutines (BLAS):
LAPACK    =-mkl
SCALAPACK =

# Gnu Scientific Library (GSL):
GSL       = -lgsl -lgslcblas
GSLHEADER =-I/usr/local/include

# Asynchronous Dynamic Load-Balancing library (ADLB)
# Note: in the default configuration, libfadlb.a, libadlb.a, and adlbf.h must be
# found in the current directory.
ADLB      =-I./ -L./ -lfadlb -ladlb 

#==============================================================================#
#                           Typical compiler flags                             #
#==============================================================================#
MACHINE = "User Defined"

ifeq ($(COMPILER),gfortran)
  # serial compiler
  FC  = gfortran
  CXX = gcc
  # Default compiler flags
  FFLAGS_OPT   =-O2
  FFLAGS_DEBUG =-O0 -Wall -Wextra -Wconversion -pedantic -fbounds-check -fimplicit-none -g -fbacktrace
  CXXFLAGS     =-O2
  OPENMP       =-fopenmp
  # suggested flags for gfortran on OS X
  ifeq ($(shell uname), Darwin)
    MACHINE = "User Defined for OS X"
    LAPACK  =-framework Accelerate
    OPENMP  =-fopenmp
    FFLAGS_OPT =-O3 -std=f2008 -pedantic
  endif
  PREPROCESSOR = -cpp -DUSE_ADLB=$(USE_ADLB)
  # tell hfbtho_basis whether basis reordering is needed
  PREPROCESSOR_BASIS = -cpp -DSORT_BASIS=$(USE_HBLAS)
endif

ifeq ($(COMPILER),ifort)
  # serial compiler
  FC  = ifort
  CXX = icc
  # Suggested flags for Intel Fortran on Linux:
  OPENMP       =-qopenmp
  FFLAGS_OPT   =-O3 -diag-disable remark -xhost -assume buffered_io -lstdc++ -stand f08
  FFLAGS_DEBUG =-O0 -diag-disable remark -debug all -check all -g -traceback -lstdc++ -stand f08
  ifeq ($(debug),1)
    CXXFLAGS =-O0 -diag-disable remark -debug all -g -traceback
  else
    CXXFLAGS =-O3 -diag-disable remark -xhost
  endif
  PREPROCESSOR = -fpp -DUSE_ADLB=$(USE_ADLB)
  # tell hfbtho_basis whether basis reordering is needed
  PREPROCESSOR_BASIS = -fpp -DSORT_BASIS=$(USE_HBLAS)
endif
 
#==============================================================================#
#                       Overide settings for specific machines                 #
#==============================================================================#
# Dogwood (UNC) - Intel
ifneq ($(shell hostname | grep dogwood),)
  MACHINE = "Dogwood"

  LAPACK    =-mkl
  SCALAPACK =

  GSL          =-L/nas/longleaf/apps/gsl/2.4/lib -lgsl -lgslcblas -cxxlib
  GSLHEADER    =-I/nas/longleaf/apps/gsl/2.4/include
endif

#==============================================================================#
#                         Compile-time preprocessing                           #
#==============================================================================#
# Consistency Check
ifeq ($(USE_ADLB),1)
  ifeq ($(USE_MPI),0)
    $(error USE_ADLB requires USE_MPI=1)
  endif
endif

# Set mpi compiler wrappers
ifeq ($(USE_MPI),1)
  FC    = mpif90
  CXX   = mpicxx
endif

# Debugging on/off
ifeq ($(debug),1)
  FFLAGS = $(FFLAGS_DEBUG)
else
  FFLAGS = $(FFLAGS_OPT)
endif

# OpenMP disabled by default but can be switched on (e.g. make openmp=1 all)
ifneq ($(openmp),1)
  OPENMP =
endif

# Unless overridden in the custom setting above, the CFLAGS should mirror the FFLAGS
ifeq ($(CXXFLAGS),)
  CXXFLAGS = $(FFLAGS)
endif

# Don't link libraries if ADLB not requested
ifneq ($(USE_ADLB),1)
  ADLB =
endif

#==============================================================================#
#                         Setup Targets and Compile                            #
#==============================================================================#
# This allows tests/Makefile to use the same settings automatically
export FC CC FFLAGS CFLAGS GSLHEADER GSL LAPACK SCALAPACK OPENMP VERBOSE DOC_DIR

### Phony Targets ---
# The targets are put together into various groups:
# all (default)      - runs nompi, mpi, beta, and hfbinfo
# mpi                - builds pnfam_parallel.x
# nompi              - builds pnfam_serial.x
# beta               - builds betadecay.x
# hfbinfo            - builds hfbinfo.x
# tests              - builds module tests
# mfam               - builds the matrix FAM solver
# clean              - removes *.o and *.mod
# cleanall           - runs clean and also removes *.x
# version            - write the current version to the source
.PHONY : all nompi mpi beta hfbinfo tests mfam clean cleanall version machine

# Main targets
all     : machine nompi mpi beta hfbinfo
mpi     : machine version pnfam_mpi.x
nompi   : machine version pnfam_nompi.x
beta    : machine version betadecay.x
hfbinfo : machine version hfbinfo.x
mfam    : machine version make_mfam

# Make the necessary modules and objects and then move to the mfam directory
machine :
	@echo "Using Machine Settings For:" $(MACHINE) "+" $(COMPILER)

make_mfam :
	@make -C mfam

tests :
	@make -C tests/modules

clean :
	@-rm -f ./*.o
	@-rm -f ./*.mod
	@-rm -f version.inc
	@cd tests/modules; make clean
	@cd mfam; make clean
	@cd doc; make clean

cleanall : clean
	@-rm -f pnfam_nompi.x pnfam_mpi.x betadecay.x hfbinfo.x
	@cd tests/modules; make cleanall
	@cd mfam; make cleanall
	@cd doc; make clean
# @cd tests; rm -f ./*/*.out ./*/*.ctr ./*/thoout.dat ./*/*.hfb ./*/*.hel ./*/*.log ./regression/test.in

# Version numbering
version :
	$(SHELL) update-version.sh

doc: pnfam_nompi.x pnfam_mpi.x
	@make -C doc

### Objects and executables ---

# Require OpenMP
ifeq ($(USE_HBLAS),1)
  h_file = hamiltonian_blas.f90
else
  h_file = hamiltonian.f90
endif
hamiltonian.o : $(h_file)
ifeq ($(VERBOSE),1)
	$(FC) -c -o $@ $< $(FFLAGS) $(OPENMP)
else
	@echo "Compiling  \"$@\"..."
	@$(FC) -c -o $@ $< $(FFLAGS) $(OPENMP)
endif
pnfam_common.o : pnfam_common.f90
ifeq ($(VERBOSE),1)
	$(FC) -c -o $@ $< $(FFLAGS) $(OPENMP)
else
	@echo "Compiling  \"$@\"..."
	@$(FC) -c -o $@ $< $(FFLAGS) $(OPENMP)
endif

# Require MPI/ADLB and OpenMP
pnfam_parallel.o : pnfam_parallel.f90
ifeq ($(VERBOSE),1)
	$(FC) -c -o $@ $< $(PREPROCESSOR) $(FFLAGS) $(OPENMP) $(ADLB)
else
	@echo "Compiling  \"$@\"..."
	@$(FC) -c -o $@ $< $(PREPROCESSOR) $(FFLAGS) $(OPENMP) $(ADLB)
endif

# Require GSL
polyfit.o : polyfit.cpp
ifeq ($(VERBOSE),1)
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(GSLHEADER)
else
	@echo "Compiling  \"$@\"..."
	@$(CXX) -c -o $@ $< $(CXXFLAGS) $(GSLHEADER)
endif

# hfbtho_basis
hfbtho_basis.o : hfbtho_basis.f90
ifeq ($(VERBOSE),1)
	$(FC) -c -o $@ $< $(FFLAGS) $(PREPROCESSOR_BASIS)
else
	@echo "Compiling  \"$@\"..."
	@$(FC) -c -o $@ $< $(FFLAGS) $(PREPROCESSOR_BASIS)
endif

# Other Fortran source
%.o : %.f90
ifeq ($(VERBOSE),1)
	$(FC) -c -o $@ $< $(FFLAGS)
else
	@echo "Compiling  \"$@\"..."
	@$(FC) -c -o $@ $< $(FFLAGS)
endif

# Other C++ source
%.o : %.cpp
ifeq ($(VERBOSE),1)
	$(CXX) -c -o $@ $< $(CXXFLAGS)
else
	@echo "Compiling  \"$@\"..."
	@$(CXX) -c -o $@ $< $(CXXFLAGS)
endif
# Version numbering
constants.o : version

# Executables
pnfam_nompi.x : constants.o logger.o blockmatrix_type.o hfbtho_basis.o extfield_type.o extfield.o interaction.o broyden_mixer.o density_set_type.o hamiltonian.o pnfam_solver.o fermi.o complex_quadrature.o polyfit.o rational_interp.o phasespace.o pnfam_common.o pnfam_serial.o pnfam_WX.o
ifeq ($(VERBOSE),1)
	$(FC) -o $@ $^ $(FFLAGS) $(OPENMP) $(LAPACK) $(GSL)
else
	@echo "Assembling  \"$@\"..."
	@$(FC) -o $@ $^ $(FFLAGS) $(OPENMP) $(LAPACK) $(GSL)
endif

pnfam_mpi.x : constants.o logger.o logger_parallel.o blockmatrix_type.o hfbtho_basis.o extfield_type.o extfield.o interaction.o broyden_mixer.o density_set_type.o hamiltonian.o pnfam_solver.o fermi.o complex_quadrature.o polyfit.o rational_interp.o phasespace.o pnfam_common.o pnfam_parallel.o
ifeq ($(VERBOSE),1)
	$(FC) -o $@ $^ $(PREPROCESSOR) $(FFLAGS) $(OPENMP) $(LAPACK) $(GSL) $(ADLB)
else
	@echo "Assembling  \"$@\"..."
	@$(FC) -o $@ $^ $(PREPROCESSOR) $(FFLAGS) $(OPENMP) $(LAPACK) $(GSL) $(ADLB)
endif

betadecay.x : constants.o blockmatrix_type.o logger.o hfbtho_basis.o complex_quadrature.o fermi.o polyfit.o rational_interp.o phasespace.o betadecay.o betadecay_main.o
ifeq ($(VERBOSE),1)
	$(FC) -o $@ $^ $(FFLAGS) $(OPENMP) $(LAPACK) $(GSL)
else
	@echo "Assembling  \"$@\"..."
	@$(FC) -o $@ $^ $(FFLAGS) $(OPENMP) $(LAPACK) $(GSL)
endif

hfbinfo.x : constants.o blockmatrix_type.o logger.o hfbtho_basis.o hfbinfo.o
ifeq ($(VERBOSE),1)
	$(FC) -o $@ $^ $(FFLAGS) $(LAPACK)
else
	@echo "Assembling  \"$@\"..."
	@$(FC) -o $@ $^ $(FFLAGS) $(LAPACK)
endif
