rootDir=./


Headers=$(wildcard *.hpp)
Path=$(PWD)


CC=g++
FLG= -O3 -std=c++17 -lm  -I$(rootDir) -Wall 

EXE = $(shell find . -type f -name '*.cpp' | sed 's/\.cpp/.run/g' | sed 's/\.\///g')

all:  $(EXE)

%.run: %.cpp makefile $(Headers) makefile 
	$(CC) -o $@ $< $(FLG) 

clean:
	rm -rf *.run
	 