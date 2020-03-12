#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.56)

AC_INIT(surrogate, 1.0.12, roconnor@health.usf.edu)
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR([gen_control_main.c])

# Checks for programs.
AC_GNU_SOURCE
AC_PROG_CC
AM_PROG_CC_C_O

# Checks for libraries.
AC_PATH_XTRA

# Checks for header files.
AC_CHECK_DECL(getline, [AC_DEFINE(HAVE_GETLINE,1,[Define to 1 if you have GNU getline])], [AC_LIBOBJ(getline)])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_SIZEOF(long)

# Checks for library functions.
AC_REPLACE_FUNCS([asprintf])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
