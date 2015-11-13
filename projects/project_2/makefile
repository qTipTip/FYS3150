CC=g++
CFLAGS=-std=c++11 -O3
SOURCES=src
TESTS=tests
BUILD=build
MAINEXE=schrodinger
TESTEXE=test

all: schrodinger test
test: test.o jacobi.o schrodinger.o
	$(CC) $(CFLAGS) -larmadillo $(BUILD)/test.o $(BUILD)/jacobi.o -o tests/$(TESTEXE)

schrodinger: jacobi.o schrodinger.o
	$(CC) $(CFLAGS) $(BUILD)/jacobi.o $(BUILD)/schrodinger.o -o build/$(MAINEXE)

jacobi.o: $(SOURCES)/jacobi.cc
	$(CC) $(CFLAGS) -c $(SOURCES)/jacobi.cc -o $(BUILD)/jacobi.o

schrodinger.o: $(SOURCES)/schrodinger.cc
	$(CC) $(CFLAGS) -c $(SOURCES)/schrodinger.cc -o $(BUILD)/schrodinger.o

test.o:
	$(CC) $(CFLAGS) -c $(TESTS)/test.cc -o $(BUILD)/test.o
clean:
	rm $(BUILD)/* $(TESTS)/test data/*

