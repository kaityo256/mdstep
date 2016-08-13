TARGET=a.out

.SUFFIXES: .c .cc .h. .o
SRC=$(shell ls *.cpp)
HED=$(shell ls *.hpp)
OBJ=$(SRC:.cpp=.o)

CC=g++
CPPFLAGS=-std=c++11 -march=native -Wall -Wextra -O3

-include makefile.opt

all: a.out

a.out: $(OBJ)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJ)

.cpp.o:
	  $(CC) $(CPPFLAGS) -c $< 

dep:
	  g++  -std=gnu++11 -MM -MG $(SRC) >makefile.depend

clean:
	rm -f $(TARGET) $(OBJ)

-include makefile.depend
