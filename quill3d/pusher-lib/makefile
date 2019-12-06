CXX := mpicxx
CXXFLAGS := -O3 -march=native -ffast-math -std=c++0x -Wall -Wextra

all: libpusher

pusher: pusher.cpp pusher.h thinparticle.h
	$(CXX) -c -o pusher.o pusher.cpp $(CXXFLAGS)

.PHONY : libpusher
libpusher: pusher.o
	ar rcs libpusher.a pusher.o

.PHONY : clean
clean:
	-rm -f pusher.o libpusher.a
