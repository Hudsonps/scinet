CC = g++-mp-6
CFLAGS = -fopenmp -m64 -std=c++11
COMPILE_COMMAND = $(CC) $(CFLAGS)
OUTPUT = my_test

all: omptest.cc
	$(COMPILE_COMMAND) -o $(OUTPUT) omptest.cc

clean:
	rm -f *.o $(OUTPUT).*
