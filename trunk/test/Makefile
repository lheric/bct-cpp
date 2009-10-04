filenames = \
	clustering_coef_bu_cpp \
	degrees_dir_cpp \
	degrees_und_cpp
sources := $(foreach f, $(filenames), $(f).cpp)
objects := $(foreach f, $(filenames), $(f).o)
oct_files := $(foreach f, $(filenames), $(f).oct)

.PHONY: all clean install uninstall

all:
	for source in $(sources); do \
		mkoctfile -lbct -lgsl $$source; \
	done

install:
	if [ ! -d /usr/local/share/bct-cpp ]; then \
		mkdir /usr/local/share/bct-cpp; \
	fi
	cp $(oct_files) /usr/local/share/bct-cpp

uninstall:
	-rm -rf /usr/local/share/bct-cpp

clean:
	-rm $(objects) $(oct_files)
