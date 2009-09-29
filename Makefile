objects = \
	degrees_dir.o \
	degrees_und.o \
	utility.o
bct: $(objects)
	libtool -static -o libbct.a $(objects)
degrees_dir.o:
	cc -c degrees_dir.cpp
degrees_und.o:
	cc -c degrees_und.cpp
utility.o:
	cc -c utility.cpp
clean:
	-rm libbct.a $(objects)
