DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MMD -MP -MF $(DEPDIR)/$*.Td
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

CXX := mpicxx
CXXFLAGS := -O3 -march=native -ffast-math -std=c++0x -Wall -Wextra
debug: CXXFLAGS := -g -std=c++0x
CXXLIBS := 
OBJECTS := advance.o deposition.o fout.o init.o low-level.o main.o mc.o plist.o maxwell.o containers.o
DEPS := $(patsubst %,$(DEPDIR)/%.d,$(basename $(OBJECTS)))

quill: $(OBJECTS) $(DEPS)
	$(MAKE) libpusher -C pusher-lib
	$(CXX) -o quill $(CXXFLAGS) $(OBJECTS) $(CXXLIBS) -Lpusher-lib -lpusher

%.o: %.cpp
%.o: %.cpp $(DEPDIR)/%.d
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -o $@ -c $<
	$(POSTCOMPILE)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

include $(wildcard $(DEPS))

debug: $(OBJECTS)
	$(CXX) -o quill $(CXXFLAGS) $(OBJECTS) $(CXXLIBS)

.PHONY : clean
clean :
	-rm -f quill $(OBJECTS) conf
	-rm -rf .d/
	$(MAKE) clean -C pusher-lib

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
	./parse.py $(FILE) > conf

.PHONY : chameleon
chameleon:
	$(MAKE) -C chameleon

all: quill chameleon

