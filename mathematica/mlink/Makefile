
# adapt path to your Mathematica installation
MATHSYSTEMDIR = /usr/applic/packages/mathematica/10.0.0/SystemFiles

# set this value to the result of evaluating $SystemID
SYS = Linux-x86-64

CADDSDIR = ${MATHSYSTEMDIR}/Links/MathLink/DeveloperKit/${SYS}/CompilerAdditions
MATHLDIR = ${MATHSYSTEMDIR}/Libraries/${SYS}

# set these with appropriate libraries for your system
EXTRALIBS = -lm -lpthread -lrt -ldl -luuid
MLLIB = ML64i4

MPREP = ${CADDSDIR}/mprep

CCOPTS = -Wall -m64 -O3

all: fermifabML

fermifabML: bitfield.o boson_map.o fermi_map.o fermifab_mlink_linux.o fermifabML.o generate_rdm.o slater_rdm.o spintrace_coul.o util.o
	${CC} ${CCOPTS} -L${CADDSDIR} -L${MATHLDIR} -l${MLLIB} ${EXTRALIBS} -o fermifabML/${SYS}/$@ $?

%.o: src/%.c
	${CC} -c -Iinclude -I${CADDSDIR} ${CCOPTS} $<

src/fermifab_mlink_linux.c: fermifab_mlink.tm
	${MPREP} -o $@ $?


clean:
	rm *.o
