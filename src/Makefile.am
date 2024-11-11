AM_CPPFLAGS = -I$(top_srcdir) -I$(includedir) -DG_LOG_DOMAIN=\"RVPM\"

lib_LTLIBRARIES = librvpm.la

librvpm_la_LDFLAGS = -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic

librvpm_la_SOURCES = \
	rvpm-private.h \
	rvpm.c \
	kernel.c \
	tree.c \
	solver.c \
	streamfunc.c \
	specfunc.c \
	common.c

include_HEADERS = \
	rvpm.h

librvpm_la_LIBADD = libsub.la
noinst_LTLIBRARIES = libsub.la
libsub_la_SOURCES = \
	rvpm-private.h \
	rvpm.c \
	kernel.c \
	tree.c \
	solver.c \
	specfunc.c \
	streamfunc.c

libsub_la_CFLAGS = -DRVPM_SINGLE_PRECISION

CLEANFILES = $(BUILT_SOURCES)
