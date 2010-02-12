CXXFLAGS = -Wall -g

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
	latmio_und.o \
	latmio_und_connected.o \
	macaque.o \
	make_motif34lib.o \
	matching_ind.o \
	matlab/compare.o \
	matlab/convert.o \
	matlab/functions.o \
	matlab/index.o \
	matlab/operators.o \
	matlab/utility.o \
	motif3funct_bin.o \
	motif3funct_wei.o \
	motif3struct_bin.o \
	motif3struct_wei.o \
	norm_avr_shortest_path_length_bu.o \
 	norm_avr_shortest_path_length_bd.o \
	norm_avr_shortest_path_length_wd.o \
	norm_avr_shortest_path_length_wu.o \
	randmio_dir.o \
	randmio_dir_connected.o \
	randmio_und.o \
	randmio_und_connected.o \
	reachdist.o \
	status.o \
	strengths_dir.o \
	strengths_und.o \
	utility.o

.PHONY: all clean install swig uninstall

all: libbct.a

libbct.a: $(objects)
	ar crs libbct.a $(objects)

swig: $(objects)
	swig -Wall -c++ -python -o bct_wrap.cpp bct.h
	g++ $(CXXFLAGS) -c -D BCT_SWIG -I/usr/include/python2.5 -include bct.h bct_wrap.cpp swig.cpp
	g++ $(CXXFLAGS) -lgsl -lgslcblas -bundle -flat_namespace -undefined suppress -o _bct.so $(objects) bct_wrap.o swig.o

$(objects): matlab/matlab.h bct.h

install: libbct.a
	if [ ! -d /usr/local/include/bct ]; then \
		mkdir /usr/local/include/bct; \
		mkdir /usr/local/include/bct/matlab; \
	fi
	cp matlab/matlab.h /usr/local/include/bct/matlab
	cp matlab/sort.h /usr/local/include/bct/matlab
	cp bct.h /usr/local/include/bct
	cp libbct.a /usr/local/lib

uninstall:
	-rm -rf /usr/local/include/bct
	-rm /usr/local/lib/libbct.a

clean:
	-rm _bct.so bct.py bct.pyc bct_wrap.cpp bct_wrap.o libbct.a swig.o $(objects)
