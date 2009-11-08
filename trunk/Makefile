objects = \
	assortativity.o \
	cat.o \
	clustering_coef_bu.o \
	clustering_coef_wu.o \
	convert.o \
	debug.o \
	degrees_dir.o \
	degrees_und.o \
	density_dir.o \
	density_und.o \
	error.o \
	fve.o \
	jdegree.o \
	macaque.o \
	matching_ind.o \
	matlab.o \
	status.o \
	strengths_dir.o \
	strengths_und.o \
	utility.o

.PHONY: all clean install uninstall

all: libbct.a

libbct.a: $(objects)
	ar crs libbct.a $(objects)

$(objects): bct.h

install: libbct.a
	cp matlab.h /usr/local/include
	cp bct.h /usr/local/include
	cp libbct.a /usr/local/lib

uninstall:
	-rm /usr/local/include/matlab.h
	-rm /usr/local/include/bct.h
	-rm /usr/local/lib/libbct.a

clean:
	-rm libbct.a $(objects)
