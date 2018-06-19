CC = mpic++
CFLAGS = -std=c++11

all:	assign2

assign2: assign2.o assign2functions.o
	$(CC) $(CFLAGS) assign2.o assign2functions.o -o assign2

assign2.o: assign2.cc
	$(CC) $(CFLAGS) -c assign2.cc

assign2functions.o: assign2functions.cc
	$(CC) $(CFLAGS) -c assign2functions.cc

wave1dorig: wave1dorig.cc
	$(CC) $(CFLAGS) wave1dorig.cc -o wave1dorig

clean:
	rm *o assign2

integratedtest:
	make ; \
	./assign2 ; \
	diff -u myresults.txt myresults_serial.txt

shapevstime_test: shapevstime_test.o assign2functions.o
	$(CC) $(FLAGS) shapevstime_test.o assign2functions.o -o shapevstime_test

unittest:
	make shapevstime_test ; \
	./shapevstime_test



