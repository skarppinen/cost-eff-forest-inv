CC:=gcc
GSL_FLAGS:=-lgsl
CBLAS_FLAGS:=-lopenblas
LIBNAME:=libpov
OPTLEVEL:=-O2 -march=native
OBJ:=random-sweep.o pov.o model-normal.o model-lognormal.o
FLAGS:=-Wall -Wextra -Wpedantic -std=c99 $(OPTLEVEL) 
.PHONY=clean install shared

all: random-sweep.o pov.o model-normal.o model-lognormal.o $(LIBNAME).so 

random-sweep.o: random-sweep.c random-sweep.h
	$(CC) $(FLAGS) -fpic -c random-sweep.c 

pov.o: pov.c pov.h
	$(CC) $(FLAGS) -fpic -c pov.c

model-normal.o: model-normal.c
	$(CC) $(FLAGS) -fpic -c model-normal.c

model-lognormal.o: model-lognormal.c
	$(CC) $(FLAGS) -fpic -c model-lognormal.c

random-sweep-prof.o: random-sweep.c random-sweep.h
	$(CC) $(FLAGS) -pg -c random-sweep.c -o random-sweep-prof.o

$(LIBNAME).so: $(OBJ)
	$(CC) $(FLAGS) -shared -o $(LIBNAME).so $(OBJ) $(GSL_FLAGS) $(CBLAS_FLAGS)
 
shared: $(LIBNAME).so

clean:
	rm *.o
	rm *.so

