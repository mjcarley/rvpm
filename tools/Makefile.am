## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)
LDADD = $(top_builddir)/src/.libs/librvpm.a 

DEPS = $(top_builddir)/src/librvpm.la

EXTRA_DIST = \
	rvpm-test.c.in \
	rvpm-solve.c.in \
	make-source \
	make-source-float

bin_PROGRAMS = \
	rvpm-test \
	rvpm-test-float \
	rvpm-solve \
	rvpm-solve-float \
	rvpm-grid \
	rvpm-util

BUILT_SOURCES = \
	rvpm-test.c \
	rvpm-test-float.c \
	rvpm-solve.c \
	rvpm-solve-float.c

CLEANFILES = $(BUILT_SOURCES)

rvpm-test.c: rvpm-test.c.in
	./make-source rvpm-test

rvpm-test-float.c: rvpm-test.c.in
	./make-source-float rvpm-test

rvpm-solve.c: rvpm-solve.c.in
	./make-source rvpm-solve

rvpm-solve-float.c: rvpm-solve.c.in
	./make-source-float rvpm-solve

