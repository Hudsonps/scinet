GSL_INC = /usr/local/include
GSL_LIB = /usr/local/lib
GSLFLAGS += "gsl-config -cflags -libs"

all:	assign5

assign5: assign5.o assign5functions.o 
	g++ -Wall -std=c++11 assign5.o assign5functions.o -L${GSL_LIB} -o assign5 -lgslcblas -lgsl

assign5.o: assign5.cc
	g++ -std=c++11 -I${GSL_INC} -c assign5.cc 

assign5functions.o: assign5functions.cc
	g++ -std=c++11 -I${GSL_INC} -c assign5functions.cc

clean:
	rm *o assign5


 

