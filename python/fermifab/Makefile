# Makefile for standalone tests

# source files
SRCFILES = src/bitfield.c src/boson_map.c src/fermi_map.c src/generate_rdm.c src/sparse.c src/tensor_op.c src/util.c
# test files
TSTFILES = test/binio.c test/test_tensor_op.c

# compiler options
CCOPTS = -Wall -O2

# set these with appropriate include paths for your system
INCLUDIRS = -Iinclude -I/usr/include/x86_64-linux-gnu

# set these with appropriate libraries for your system;
# the CBLAS and LAPACKE libraries can be installed via `sudo apt install libblas-dev liblapacke-dev` (on Ubuntu Linux) or similar
LIBRARIES = -lm -lblas -llapacke


all: proj_test

proj_test: ${SRCFILES} ${TSTFILES} test/run_tests.c
	${CC} ${CCOPTS} ${INCLUDIRS} -g -o test/run_tests $? ${LIBRARIES}
