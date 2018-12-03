# C compiler and flags
CPPFLAGS = -I/opt/local/include
LDFLAGS = -L/opt/local/lib -lgsl -lgslcblas -lfftw3f_omp -lfftw3f -lm
CC = /usr/local/bin/gcc -fopenmp -Ofast
#CC = gcc -fopenmp

UNAME := $(shell uname)

# object files
OBJ_FILES = Compute_21cmPS

#########################################################################

ifeq ($(UNAME), Linux)
Compute_21cmPS: Compute_21cmPS.c \
	${COSMO_FILES} \
	
	${CC} ${CPPFLAGS} -fPIC -c Compute_21cmPS.c ${LDFLAGS}
	
	${CC} -shared -Wl,-soname,Compute_21cmPS.so.1 -o Compute_21cmPS.so Compute_21cmPS.o ${LDFLAGS}
endif
ifeq ($(UNAME), Darwin)
Compute_21cmPS: Compute_21cmPS.c \
        ${COSMO_FILES} \

	${CC} ${CPPFLAGS} -o Compute_21cmPS.so -shared Compute_21cmPS.c ${LDFLAGS}
endif

clean:
	rm *.o *.so
