rootDir=./


Headers=$(wildcard *.hpp)
Path=$(PWD)


CC=g++
FLG= -O3 -std=c++17 -lm  -I$(rootDir) -Wall 


all:  1D-example.run 3D-example.run

1D-example.run: 1D-example.cpp  $(Headers) makefile 
	$(CC) -o $@ $< $(FLG) 

3D-example.run: 3D-example.cpp  $(Headers) makefile 
	$(CC) -o $@ $< $(FLG) 


clean:
	rm -rf *.run
	 