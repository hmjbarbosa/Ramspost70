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
LIBUTILS=$(UTILS_LIB)/libutils-g95-$(UTILS_VERSION)-ramspost.a

# Machine-dependent options.
#---------------------  LINUX Gnu gcc/g95  -----------------------
CMACH=PC_LINUX1
F_COMP=g95
F_OPTS=-O3 -fno-second-underscore -fmultiple-save
C_COMP=gcc
C_OPTS=-O3 -DUNDERSCORE -DLITTLE
LOADER=g95
#LOADER_OPTS=-O3 -fno-second-underscore
LOADER_OPTS=-O3
LIBS=
#-----------------------------------------------------------------

# For IBM,HP,SGI,ALPHA use these:
#ARCHIVE=ar rs
# For SUN,CONVEX, LINUX
ARCHIVE=ar r
