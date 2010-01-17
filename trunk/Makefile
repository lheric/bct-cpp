objects = \
	assortativity.o \
	betweenness_bin.o \
	betweenness_wei.o \
	breadth.o \
	breadthdist.o \
	cat.o \
	charpath.o \
	clustering_coef_bd.o \
	clustering_coef_bu.o \
	clustering_coef_wd.o \
	clustering_coef_wu.o \
	convert.o \
	cycprob.o \
	debug.o \
	degrees_dir.o \
	degrees_und.o \
	density_dir.o \
	density_und.o \
	distance_bin.o \
	distance_wei.o \
	efficiency.o \
	erange.o \
	findpaths.o \
	findwalks.o \
	fve.o \
	jdegree.o \
	latmio_dir.o \
	latmio_dir_connected.o \
	latmio_und_connected.o \
	latmio_und.o \
	randmio_und_connected.o \
	macaque.o \
	make_motif34lib.o \
	matching_ind.o \
	matlab/compare.o \
	matlab/convert.o \
	matlab/functions.o \
	matlab/index.o \
	matlab/operators.o \
	norm_avr_shortest_path_length_bu.o \
	reachdist.o \
	status.o \
	strengths_dir.o \
	strengths_und.o \
	utility.o

.PHONY: all clean install uninstall

all: libbct.a

libbct.a: $(objects)
	ar crs libbct.a $(objects)

$(objects): matlab/matlab.h bct.h

install: libbct.a
	if [ ! -d /usr/local/include/bct ]; then \
		mkdir /usr/local/include/bct; \
		mkdir /usr/local/include/bct/matlab; \
	fi
	cp matlab/matlab.h /usr/local/include/bct/matlab
	cp matlab/quicksort.h /usr/local/include/bct/matlab
	cp bct.h /usr/local/include/bct
	cp libbct.a /usr/local/lib

uninstall:
	-rm -rf /usr/local/include/bct
	-rm /usr/local/lib/libbct.a

clean:
	-rm libbct.a $(objects)
