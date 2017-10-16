CPP := g++
CPPOPTS := -O3 -std=c++11 -shared -fPIC
# path for R.h and Rinternals.h
#RPATH := /usr/include/R
RPATH := /usr/share/R/include
# directory of libchameleon, and linker options
LOPTS := -L. -Wl,-rpath=.

# see for compiling and use of shared objects:
# http://www.cprogramming.com/tutorial/shared-libraries-linux-gcc.html

libchameleon.so: chameleon-core.cpp
	$(CPP) -o $@ $(CPPOPTS) chameleon-core.cpp

libchameleon-R.so: chameleon-R.cpp libchameleon.so
	$(CPP) -o $@ $(CPPOPTS) -I$(RPATH) $(LOPTS) chameleon-R.cpp -lchameleon

clean:
	rm -f libchameleon.so libchameleon-R.so

example-py: libchameleon.so
	python chameleon-example.py

example-R: libchameleon-R.so
	R --slave -f chameleon-example.R
