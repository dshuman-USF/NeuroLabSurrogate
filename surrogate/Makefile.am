#                                               -*- Makefile -*-

bin_PROGRAMS = gen_control unpin gammafit_search surrogate	\
edt_surrogate explain miura eval checksurr repat

eval_SOURCES = eval.c edt.c edt.h

repat_SOURCES = repat.c edt.c edt.h

gen_control_SOURCES = gen_control_main.c gen_control.c read_array.c	\
gen_control.h read_array.h fpu_check.c fpu_check.h

unpin_SOURCES = unpin.c unpin_main.c unpin.h fpu_check.c fpu_check.h	\
unpin_nopin.c

gammafit_search_SOURCES = gammafit_search.c gammafit_search_main.c	\
gammafit_search.h read_array.c read_array.h fpu_check.c fpu_check.h	\
gammafit_search_nopin.c util.h util.c

explain_SOURCES = explain.c unpin.c unpin.h fpu_check.c fpu_check.h	\
unpin_nopin.c

surrogate_SOURCES = surrogate.c read_array.c read_array.h	\
gen_control.c gen_control.h unpin.c unpin.h gammafit_search.c	\
gammafit_search.h fpu_check.c fpu_check.h unpin_nopin.c		\
gammafit_search_nopin.c util.c util.h

edt_surrogate_SOURCES = edt_surrogate.c edt.c edt.h gen_control.c	\
gen_control.h unpin.c unpin.h gammafit_search.c gammafit_search.h	\
fpu_check.c fpu_check.h unpin_nopin.c gammafit_search_nopin.c util.c	\
util.h tools.c tools.h

miura_SOURCES = miura.c edt.c edt.h tools.c tools.h

checksurr_SOURCES = checksurr.c edt.c edt.h tools.c tools.h

bin_SCRIPTS = dogamma

EXTRA_DIST = dogamma spikes-1-0.1.gdf spikes-1-0.2.gdf			\
spikes-1-0.4.gdf spikes-16-0.1.gdf spikes-16-0.2.gdf spikes-16-0.4.gdf	\
spikes-2-0.1.gdf spikes-2-0.2.gdf spikes-2-0.4.gdf spikes-30-0.1.gdf	\
spikes-30-0.2.gdf spikes-30-0.4.gdf spikes-4-0.1.gdf spikes-4-0.2.gdf	\
spikes-4-0.4.gdf spikes-8-0.1.gdf spikes-8-0.2.gdf spikes-8-0.4.gdf	\
dogamma.out all103 gen_control_103_26_0.out debian

LDADD =$(LIBOBJS) -lgsl -lgslcblas -lJudy -lm
AM_CFLAGS = -Wall --std=c99

explain_CFLAGS = -Wall $(X_CFLAGS)
explain_LDADD = $(LIBOBJS) $(X_LIBS) -lX11 -lJudy -lm -lgsl -lgslcblas

dist_doc_DATA = README
