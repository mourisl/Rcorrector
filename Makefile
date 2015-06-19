CXX = g++
CXXFLAGS= -Wall -O3 -std=c++0x
LINKFLAGS = -I . -lpthread  
DEBUG=
OBJECTS = KmerCode.o ErrorCorrection.o
all: main.o $(OBJECTS)
	$(CXX) -o rcorrector $(CXXFLAGS) $(OBJECTS) main.o $(LINKFLAGS) 

main.o: main.cpp utils.h Reads.h Store.h
KmerCode.o: KmerCode.cpp KmerCode.h 
ErrorCorrection.o: ErrorCorrection.cpp ErrorCorrection.h

clean:
	rm *.o *.gch rcorrector
