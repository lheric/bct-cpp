CXXFLAGS = -Wall -g -arch i386

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
	connectivity_length.o \
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
	find_motif34.o \
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
	makeevenCIJ.o \
	makefractalCIJ.o \
	makelatticeCIJ.o \
	makerandCIJ_bd.o \
	makerandCIJ_bu.o \
	makerandCIJ_wd.o \
	makerandCIJ_wu.o \
	makerandCIJdegreesfixed.o \
	makeringlatticeCIJ.o \
	maketoeplitzCIJ.o \
	matching_ind.o \
	matlab/matlab.o \
	matlab/matlab_double.o \
	matlab/matlab_float.o \
	matlab/matlab_long_double.o \
	modularity_louvain.o \
	modularity_newman.o \
	module_degree_zscore.o \
	motif3funct_bin.o \
	motif3funct_wei.o \
	motif3struct_bin.o \
	motif3struct_wei.o \
	motif4funct_bin.o \
	motif4funct_wei.o \
	motif4struct_bin.o \
	motif4struct_wei.o \
	normalized_path_length.o \
	participation_coef.o \
	randmio_dir.o \
	randmio_dir_connected.o \
	randmio_und.o \
	randmio_und_connected.o \
	reachdist.o \
	status.o \
	strengths_dir.o \
	strengths_und.o \
	threshold_absolute.o \
	threshold_proportional.o \
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
	cp matlab/macros.h /usr/local/include/bct/matlab
	cp matlab/matlab.h /usr/local/include/bct/matlab
	cp matlab/matlab_double.h /usr/local/include/bct/matlab
	cp matlab/matlab_float.h /usr/local/include/bct/matlab
	cp matlab/matlab_long_double.h /usr/local/include/bct/matlab
	cp matlab/sort.h /usr/local/include/bct/matlab
	cp bct.h /usr/local/include/bct
	cp libbct.a /usr/local/lib

uninstall:
	-rm -rf /usr/local/include/bct
	-rm /usr/local/lib/libbct.a

clean:
	-rm _bct.so bct.py bct.pyc bct_wrap.cpp bct_wrap.o libbct.a swig.o $(objects)
