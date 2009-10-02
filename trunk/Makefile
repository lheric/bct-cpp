objects = \
	cat.o \
	clustering_coef_bu.o \
	degrees_dir.o \
	degrees_und.o \
	fve.o \
	macaque.o \
	utility.o
bct: $(objects)
	ar crs libbct.a $(objects)
cat.o:
	cc -c cat.cpp
clustering_coef_bu.o:
	cc -c clustering_coef_bu.cpp
degrees_dir.o:
	cc -c degrees_dir.cpp
degrees_und.o:
	cc -c degrees_und.cpp
fve.o:
	cc -c fve.cpp
macaque.o:
	cc -c macaque.cpp
utility.o:
	cc -c utility.cpp
clean:
	-rm libbct.a $(objects)
