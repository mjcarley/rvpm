AC_INIT([rvpm],1.0.0)
AC_CONFIG_SRCDIR([src/rvpm.h])

RVPM_MAJOR_VERSION=1
RVPM_MINOR_VERSION=0
RVPM_MICRO_VERSION=0
RVPM_INTERFACE_AGE=0
RVPM_BINARY_AGE=0
RVPM_VERSION=$RVPM_MAJOR_VERSION.$RVPM_MINOR_VERSION.$RVPM_MICRO_VERSION

AC_SUBST(RVPM_MAJOR_VERSION)
AC_SUBST(RVPM_MINOR_VERSION)
AC_SUBST(RVPM_MICRO_VERSION)
AC_SUBST(RVPM_VERSION)

# libtool versioning
LT_RELEASE=$RVPM_MAJOR_VERSION.$RVPM_MINOR_VERSION
LT_CURRENT=`expr $RVPM_MICRO_VERSION - $RVPM_INTERFACE_AGE`
LT_REVISION=$RVPM_INTERFACE_AGE
LT_AGE=`expr $RVPM_BINARY_AGE - $RVPM_INTERFACE_AGE`
AC_SUBST(LT_RELEASE)
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)

# For automake.
VERSION=$RVPM_VERSION
PACKAGE=rvpm

AC_SUBST(PACKAGE)
AC_SUBST(VERSION)

AM_INIT_AUTOMAKE(subdir-objects)

##$PACKAGE, $VERSION, no-define)

# Specify a configuration file
AM_CONFIG_HEADER(config.h)

AC_DEFINE_UNQUOTED(RVPM_MAJOR_VERSION, $RVPM_MAJOR_VERSION, [Major version])
AC_DEFINE_UNQUOTED(RVPM_MINOR_VERSION, $RVPM_MINOR_VERSION, [Minor version])
AC_DEFINE_UNQUOTED(RVPM_MICRO_VERSION, $RVPM_MICRO_VERSION, [Micro version])
AC_DEFINE_UNQUOTED(RVPM_INTERFACE_AGE, $RVPM_INTERFACE_AGE, [Interface age])
AC_DEFINE_UNQUOTED(RVPM_BINARY_AGE, $RVPM_BINARY_AGE, [Binary age])

dnl Initialize libtool
AM_PROG_LIBTOOL

dnl Initialize maintainer mode
#AM_MAINTAINER_MODE

AC_PROG_CC
AM_PROG_CC_C_O

if test x$GCC = xyes ; then
  CFLAGS="$CFLAGS -Wall -Werror-implicit-function-declaration -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations"
fi

AC_PROG_AWK
AC_ISC_POSIX
AC_HEADER_STDC


PKG_CHECK_MODULES([GLIB], [glib-2.0], [],
[
  echo \#################################################################
  echo \# 
  echo \# GLIB-2.0 not installed. You can install it using your 
  echo \# distribution's package manager or from:
  echo \#
  echo \# https://gitlab.gnome.org/GNOME/glib/
  echo \#
  echo \#################################################################
  exit
])

PKG_CHECK_MODULES([BLASWRAP], [blaswrap], [],
[
  echo \#################################################################
  echo \#
  echo \# BLAS wrapper library not installed or not properly installed.
  echo \# You can get it from
  echo \#
  echo \# https://github.com/mjcarley/blaswrap
  echo \#
  echo \#################################################################
  exit
])

PKG_CHECK_MODULES([GRBF], [grbf], [],
[
  echo \#################################################################
  echo \# 
  echo \# GRBF \(Gaussian Radial Basis Function library\) not
  echo	\# installed. You can get it from
  echo \#
  echo \# https://github.com/mjcarley/grbf
  echo \#
  echo \#################################################################
  exit
])

PKG_CHECK_MODULES([WBFMM], [wbfmm], [],
[
  echo \#################################################################
  echo \# 
  echo \# WBFMM \(Wide Band Fast Multipole Method\) library not
  echo \# installed. You can get it from
  echo \#
  echo \# https://github.com/mjcarley/wbfmm
  echo \#
  echo \#################################################################
  exit
])

PKG_CHECK_MODULES([GQR], [gqr], [],
[
  echo \#################################################################
  echo \# 
  echo \# GQR \(Gaussian type Quadrature Rules\) library not
  echo \# installed. You can get it from
  echo \#
  echo \# https://github.com/mjcarley/gqr
  echo \#
  echo \#################################################################
  exit
])

AX_GCC_X86_CPU_SUPPORTS(avx, hasavx=yes) 
AX_GCC_X86_CPU_SUPPORTS(avx2, hasavx2=yes) 
AX_GCC_X86_CPU_SUPPORTS(fma, hasfma=yes) 

if test x$hasavx = xyes; then
   SIMD_FLAGS="$SIMD_FLAGS -mavx -DWBFMM_USE_AVX"
fi
if test x$hasfma = xyes; then
   SIMD_FLAGS="$SIMD_FLAGS -mfma -DWBFMM_USE_AVX"
fi

AC_CHECK_LIB(m, cos)
AC_CHECK_LIB(matheval, evaluator_evaluate)
AC_CONFIG_MACRO_DIR([m4])

AC_OPENMP
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN
AC_F77_WRAPPERS

#DX_DOXYGEN_FEATURE(ON)
##DX_DOT_FEATURE(ON)
#DX_HTML_FEATURE(ON)
##DX_CHM_FEATURE(OFF)
##DX_CHI_FEATURE(OFF)
##DX_MAN_FEATURE(OFF)
##DX_RTF_FEATURE(OFF)
##DX_XML_FEATURE(OFF)
##DX_PDF_FEATURE(OFF)
##DX_PS_FEATURE(OFF)
#DX_INIT_DOXYGEN(rvpm, doc/rvpm.dxy)

CPPFLAGS="$glib_cflags $blaswrap_cflags"
CFLAGS="$CFLAGS $OPENMP_CFLAGS $SIMD_FLAGS $GLIB_CFLAGS $GRBF_CFLAGS $WBFMM_CFLAGS $GQR_CFLAGS $blaswrap_cflags"
LIBS="$GLIB_LIBS $GRBF_LIBS $WBFMM_LIBS $GQR_LIBS $blaswrap_libs $LIBS"

AC_SUBST(CPPFLAGS)
AC_SUBST(CFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(XFOIL)

LIBS="$glib_libs $LIBS $FLIBS"

AC_CONFIG_FILES([
Makefile
rvpm.pc
src/Makefile
tools/Makefile])

AC_OUTPUT

