CXX = g++
CXXFLAGS= -Wall -O3 -std=c++0x
LINKFLAGS = -I . -lpthread -lz 
DEBUG=
jellyfishVersion = `jellyfish --version | cut -f2 -d' ' | cut -f1 -d'.'`
OBJECTS = KmerCode.o ErrorCorrection.o
all: main.o $(OBJECTS)
	$(CXX) -o rcorrector $(CXXFLAGS) $(OBJECTS) main.o $(LINKFLAGS)
	if [ $(jellyfishVersion) -ne 2 ] ; \
	then \
		if [ ! -f ./jellyfish/bin/jellyfish ] ; \
		then \
			wget -O jellyfish.tar.gz https://github.com/gmarcais/Jellyfish/releases/download/v2.2.4/jellyfish-2.2.4.tar.gz ;\
			tar -xzf jellyfish.tar.gz ; mv jellyfish-2.2.4 jellyfish ; rm jellyfish.tar.gz ; cd jellyfish && ./configure && make ;\
		fi;\
	fi ;

main.o: main.cpp utils.h Reads.h Store.h
KmerCode.o: KmerCode.cpp KmerCode.h 
ErrorCorrection.o: ErrorCorrection.cpp ErrorCorrection.h

clean:
	rm *.o *.gch rcorrector
