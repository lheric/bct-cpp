objects = \
	binary_matrix.o \
	binary_vector.o \
	degrees_dir.o \
	degrees_und.o \
	matrix_printf.o \
	vector_printf.o
bct: $(objects)
	libtool -static -o libbct.a $(objects)
binary_matrix.o:
	cc -c binary_matrix.cpp
binary_vector.o:
	cc -c binary_vector.cpp
degrees_dir.o:
	cc -c degrees_dir.cpp
degrees_und.o:
	cc -c degrees_und.cpp
matrix_printf.o:
	cc -c matrix_printf.cpp
vector_printf.o:
	cc -c vector_printf.cpp
clean:
	-rm libbct.a $(objects)
