CXX = g++
CXXFLAGS= -Wall -lpthread -g
DEBUG=
OBJECTS = KmerCode.o ErrorCorrection.o
all: main.o $(OBJECTS)
	$(CXX) -o rcorrector $(CXXFLAGS) $(OBJECTS) main.o 
main.o: main.cpp utils.h Reads.h Store.h
#bloom_filter.o: bloom_filter.h 
#Store.o: Store.h 
#BitTable.o: BitTable.cpp BitTable.h
KmerCode.o: KmerCode.cpp KmerCode.h 
ErrorCorrection.o: ErrorCorrection.cpp ErrorCorrection.h

clean:
	rm *.o *.gch
