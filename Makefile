CC=g++
CFLAGS=-Wall -Werror -O3 -std=c++17
CINCLUDE=-I/usr/include/python3.6
CLIBS=-lpython3.6m
BUILD_DIR=build

all: optman

optman: $(BUILD_DIR)/main.o
	cd $(BUILD_DIR); $(CC) main.o -o optman.out $(CFLAGS) $(CLIBS); cd ../

$(BUILD_DIR)/main.o: main.cpp
	$(CC) -c main.cpp -o $(BUILD_DIR)/main.o $(CFLAGS) $(CINCLUDE)

clean:
	cd $(BUILD_DIR); rm -f *.out *.o; cd ../