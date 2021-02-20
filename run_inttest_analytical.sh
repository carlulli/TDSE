
IDIR="./include"
MDIR="./modules"
TDIR="./test"
# SDIR="./scripts"

# INCKISS

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS="-g -Wall -Werror=implicit-function-declaration -fPIC"

# the build target executable:
TARGET="inttest_analytical"

# inclde directory
INCLUDE="-I include -I kissff"

# lib
# LIBKISSFFT = kissfft
# LIBS = -L $(LIBKISSFFT)

# include modules
MODULES="${MDIR}/wavefunction.c ${MDIR}/integrator.c ${MDIR}/geometry.c ${MDIR}/linearalgebra.c ${MDIR}/hamiltonian.c ${MDIR}/conjugategradient.c"

# cd ..
gcc ${CFLAGS} ${INCLUDE} ${MODULES} ${TDIR}/${TARGET}.c -Lkissfft -lkissfft-double -o ${TARGET}.exe -lm
