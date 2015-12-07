CPP := g++
CPPOPTS := -O3 -std=c++11 -shared -fPIC

chameleon-core.so: chameleon-core.cpp
	$(CPP) -o $@ $(CPPOPTS) chameleon-core.cpp

clean:
	rm chameleon-core.so

example: chameleon-core.so
	python chameleon-example.py
