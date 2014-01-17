CFLAGS=-std=c++0x -g $(shell pkg-config --cflags liblcb-experimental) 
LDFLAGS=$(shell pkg-config --libs liblcb-experimental) -lboost_program_options -lboost_system -lboost_filesystem 

OBJECTS = main.o \
  Support.o \
  Protein.o \
  Normal.o  \
  VonMises3D.o \
  Component.o \
  Mixture.o

all: main 

main: $(OBJECTS)
	g++ $(OBJECTS) -o $@ $(LDFLAGS) 

main.o: main.cpp Support.h Header.h
	g++ -c $(CFLAGS) $< -o $@

Support.o: Support.cpp Support.h Header.h
	g++ -c $(CFLAGS) $< -o $@

Protein.o: Protein.cpp Protein.h Support.h 
	g++ -c $(CFLAGS) $< -o $@

Normal.o: Normal.cpp Normal.h 
	g++ -c $(CFLAGS) $< -o $@

VonMises3D.o: VonMises3D.cpp VonMises3D.h 
	g++ -c $(CFLAGS) $< -o $@

Component.o: Component.cpp Component.h Header.h 
	g++ -c $(CFLAGS) $< -o $@

Mixture.o: Mixture.cpp Mixture.h Header.h 
	g++ -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o *~ main 

