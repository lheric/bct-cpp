filenames = \
	assortativity_dir_cpp \
	assortativity_und_cpp \
	clustering_coef_bd_cpp \
	clustering_coef_bu_cpp \
	clustering_coef_wd_cpp \
	clustering_coef_wu_cpp \
	degrees_dir_cpp \
	degrees_und_cpp \
	density_dir_cpp \
	density_und_cpp \
	jdegree_cpp \
	jdegree_id_cpp \
	jdegree_od_cpp \
	jdegree_bl_cpp \
	matching_ind_cpp \
	strengths_und_cpp \
	strengths_dir_cpp \

objects := $(foreach f, $(filenames), $(f).o)
oct_files := $(foreach f, $(filenames), $(f).oct)
m_files = \
	bct_test.m \
	bct_test_all.m

.PHONY: all clean install uninstall

all: $(oct_files)

.SECONDEXPANSION:

$(oct_files): ../libbct.a bct_test.h bct_test.cpp $$(basename $$@).cpp
	@if [ $(subst _cpp,,../$(basename $@).cpp) -nt $@ ]; then \
		echo Making $@...; \
		mkoctfile -lbct -lgsl -lgslcblas $(basename $@).cpp; else \
		echo "make: \`$@' is up to date."; \
		touch $@; \
	fi

install: $(oct_files)
	if [ ! -d /usr/local/share/bct ]; then \
		mkdir /usr/local/share/bct; \
	fi
	cp $(oct_files) $(m_files) /usr/local/share/bct

uninstall:
	-rm -rf /usr/local/share/bct

clean:
	-rm $(objects) $(oct_files)
