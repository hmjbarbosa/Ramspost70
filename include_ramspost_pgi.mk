#Makefile include include.mk
############################## Change Log ##################################
# 1.0.0.1
#
# 000823 MJB include.mk-mrc ##
#            New - defines all make environment varaibles and is included
#            in all make files. ##
#
############################################################################

# RAMS root directory.
RAMSPOST=.
UTILS_VERSION=2.0

# Source directories.
UTILS_LIB=$(RAMSPOST)/LIB
UTILS_INCS=$(RAMSPOST)/include

# MRC libraries.
LIBUTILS=$(UTILS_LIB)/libutils-pgi-$(UTILS_VERSION)-ramspost.a

# Machine-dependent options.
#-----------------  LINUX Portland Group pgf77/gcc ---------------
CMACH=PC_LINUX1
F_COMP=pgf90
F_OPTS=-Mvect=cachesize:524288,sse -Mcache_align -Munroll -Mnoframe -O2 -pc 64
C_COMP=gcc
C_OPTS=-O3 -DUNDERSCORE -DLITTLE
LOADER=pgf90
LOADER_OPTS=-v -Wl,-static
LIBS=
#-----------------------------------------------------------------

# For IBM,HP,SGI,ALPHA use these:
#ARCHIVE=ar rs
# For SUN,CONVEX, LINUX
ARCHIVE=ar r
