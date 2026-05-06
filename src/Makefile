FC=mpif90

LDFLAGS= -Wl,--start-group ${MKL_LIB}/libmkl_blas95_lp64.a ${MKL_LIB}/libmkl_lapack95_lp64.a \
        ${MKL_LIB}/libmkl_scalapack_lp64.a ${MKL_LIB}/libmkl_intel_lp64.a ${MKL_LIB}/libmkl_core.a \
        ${MKL_LIB}/libmkl_sequential.a  ${MKL_LIB}/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group \
        -L${HDF5_LIB} -lhdf5 -lpthread -lm -ldl 
FFLAGS= -O3 -g -heap-arrays -warn noalignments -I${MKL_INC}/intel64/lp64 -I${MKL_INC} -I${HDF5_INC} -I$ZLIB_INC

SOURCES=${wildcard *.f90}
OBJECTS=${SOURCES:.f90=.o}

%.o:%.f90
	${FC} -c ${FFLAGS} $^


midas: ${OBJECTS}
	${FC} -o midas ${OBJECTS} ${FFLAGS} ${LDFLAGS} 

clean:
	@-rm *.o *.mod midas

modules:
	${FC} -c ${FFLAGS} Initialization.f90 Common_variables.f90 stack_structure.f90 nrtype.f90 DiverseUtil.f90 precmod.f90 stringmod.f90
