HEADER_PATH = ./headers
CC          = g++
CFLAGS      = --std=c++11 -Wall -I$(HEADER_PATH) -c
vpath %.cpp source
vpath %.hpp headers

all: annealing relaxation

annealing: annealing.o basicFunctions.o graphFunctions.o ran2.o
	$(CC) $< basicFunctions.o graphFunctions.o ran2.o -o $@

relaxation: relaxation.o basicFunctions.o graphFunctions.o ran2.o
	$(CC) $< basicFunctions.o graphFunctions.o ran2.o -o $@

annealing.o: annealing.cpp 
	$(CC) $(CFLAGS) $<

relaxation.o: relaxation.cpp 
	$(CC) $(CFLAGS) $<

basicFunctions.o: basicFunctions.cpp basicFunctions.hpp
	$(CC) $(CFLAGS) $<

graphFunctions.o: graphFunctions.cpp graphFunctions.hpp
	$(CC) $(CFLAGS) $<

ran2.o: ran2.cpp ran2.hpp
	$(CC) $(CFLAGS) $<
	

clean:
	rm *.o annealing relaxation