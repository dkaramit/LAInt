rootDir=./


Headers=$(wildcard *.hpp)
Path=$(PWD)


CC=g++
FLG= -O3 -std=c++17 -lm  -I$(rootDir) -Wall 


all:  example.run

example.run: example.cpp  $(Headers) makefile 
	$(CC) -o $@ $< $(FLG) 


clean:
	rm -rf *.run
	 