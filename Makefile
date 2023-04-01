CC=mpic++
CFLAGS=-Wall -std=c++17 -O3 -funroll-loops

UNAME := $(shell uname)
INCLUDES=
LDFLAGS=
ifeq ($(UNAME), Darwin)
	INCLUDES=-I/opt/homebrew/opt/libomp/include
	LDFLAGS=-L/opt/homebrew/opt/libomp/lib -lomp
endif

main:
	$(CC) $(CFLAGS) $(INCLUDES) code/main.cpp $(LDFLAGS) -o a3

nosort:
	$(CC) $(CFLAGS) $(INCLUDES) code/main_nosort.cpp $(LDFLAGS) -o a3_nosort

clean:
	rm a3
