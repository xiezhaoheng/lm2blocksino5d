bin = ../bin/lm2blocksino_5d
inc = ../../include

CXX = g++
CXXFLAGS = -std=c++14

$(bin) : main.o
	$(CXX) $(CXXFLAGS) -o $(bin) main.o

main.o : main.cpp

.PHONY : clean

clean :
	-rm -v $(bin) main.o
