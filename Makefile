#!/bin/bash

#Directory containing cpp and h files from ESU
DIR_ESU= src/

CC=g++ -std=c++0x
EXEC= exec_esu

#Enabling OpenMP
OMP= -fopenmp

SRC_ESU=\
			$(DIR_ESU)main.cpp				\
			$(DIR_ESU)Esu.cpp					\
			$(DIR_ESU)Isomorphism.cpp	\
			$(DIR_ESU)Esu_graph.cpp

OBJ_ESU = ${SRC_ESU:.cpp=.o}


all: $(EXEC)

exec_esu : $(OBJ_ESU)
	$(CC) $(OMP) -o $@ $^

##############################################

%.o: %.cpp
	$(CC) $(OMP) -c -o $@ $+

clean:
	rm $(DIR_ESU)*.o
	rm $(EXEC)
