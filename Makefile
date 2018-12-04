CC=g++
CFLAGS=-Wall -Werror -O3 -std=c++17
CINCLUDE=-I/usr/include/python3.6
CLIBS=-lpython3.6m

all: optman

optman: main.o
	$(CC) main.o -o optman.out $(CFLAGS) $(CLIBS)

main.o: main.cpp
	$(CC) -c main.cpp $(CFLAGS) $(CINCLUDE)

clean:
	rm -f *.out *.o