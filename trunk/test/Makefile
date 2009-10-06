filenames = \
	clustering_coef_bu_cpp \
	degrees_dir_cpp \
	degrees_und_cpp \
	density_dir_cpp \
	density_und_cpp \
	strengths_und_cpp
objects := $(foreach f, $(filenames), $(f).o)
oct_files := $(foreach f, $(filenames), $(f).oct)
m_files = \
	bct_test.m \
	bct_test_all.m

.PHONY: all clean install uninstall

all: $(oct_files)

.SECONDEXPANSION:

$(oct_files): bct_test.h bct_test.cpp $$(basename $$@).cpp
	mkoctfile -lbct -lgsl $(basename $@).cpp

install: $(oct_files)
	if [ ! -d /usr/local/share/bct ]; then \
		mkdir /usr/local/share/bct; \
	fi
	cp $(oct_files) $(m_files) /usr/local/share/bct

uninstall:
	-rm -rf /usr/local/share/bct-cpp

clean:
	-rm $(objects) $(oct_files)
