CPP := g++
CPPOPTS := -O3 -std=c++11 -shared -fPIC
# directory of libchameleon, and linker options
LOPTS := -L. -Wl,-rpath=.

# see for compiling and use of shared objects:
# http://www.cprogramming.com/tutorial/shared-libraries-linux-gcc.html

libchameleon.so: chameleon-core.cpp
	$(CPP) -o $@ $(CPPOPTS) chameleon-core.cpp

clean:
	rm -f libchameleon.so 

example-py: libchameleon.so
	python chameleon-example.py
