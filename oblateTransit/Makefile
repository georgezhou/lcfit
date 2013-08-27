OBJECT = oblateness.o main.o
PROJECT = oblateness
SOURCES= oblateness.cc oblateness.h
INTERFACES = oblateness.i
WRAPPERS = $(INTERFACES:.i=_wrap.cxx)
PROXIES = $(INTERFACES:.i=.py)

.PHONY: all
all: $(WRAPPERS) $(SOURCES)
	./setup.py build_ext -i

%_wrap.cxx: %.i %.h ./numpy.i
	swig -c++ -python $<


.PHONY: cc
cc: $(PROJECT)

$(PROJECT) :$(OBJECT)
	$(CXX) -pg -o $@ $^

.PHONY: clean
clean:
	$(RM) $(PROJECT)
	$(RM) *.so *.pyc *_wrap.h
	$(RM) -r *.dSYM
	$(RM) -r build
	$(RM) $(WRAPPERS)
	$(RM) $(PROXIES)
	$(RM) *.o
	$(RM) .depend
	$(RM) *.out
	$(RM) *~
