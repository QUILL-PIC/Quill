CXX := g++
CXXFLAGS := -O3 -march=native -ffast-math -std=c++0x
debug: CXXFLAGS := -g -std=c++0x
CXXLIBS := -pthread -lnuma
OBJECTS := advance.o deposition.o fout.o init.o listener.o low-level.o main.o mc.o plist.o

quill: $(OBJECTS)
	$(CXX) -o quill $(CXXFLAGS) $(OBJECTS) $(CXXLIBS)
advance.o : main.h
deposition.o : main.h
fout.o : main.h
init.o : main.h
listener.o : main.h
low-level.o : main.h
main.o : main.h
mc.o : main.h
plist.o : main.h

debug: $(OBJECTS)
	$(CXX) -o quill $(CXXFLAGS) $(OBJECTS) $(CXXLIBS)

parallel: parallel.cpp
	$(CXX) -o parallel $(CXXFLAGS) parallel.cpp -pthread

.PHONY : clean
clean:
	-rm -f quill $(OBJECTS) conf parallel

.PHONY : results_clean
results_clean:
	-rm -I results/*

.PHONY : asy_clean
asy_clean:
	-rm -I asy/*.eps asy/movie/*.png asy/movie/*.gif

.PHONY : py_clean
py_clean:
	-rm -I python/*.png python/movie.gif python/*.pyc

# config for debugging in gdb by 'run < conf'
# example: make conf FILE='.beam-lp-colliding'
FILE := ""
.PHONY : conf
conf:
	./parse.sh $(FILE) > conf
