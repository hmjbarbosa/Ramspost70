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
LIBUTILS=$(UTILS_LIB)/libutils-af95-$(UTILS_VERSION)-ramspost.a

# Machine-dependent options.
#---------------------  LINUX Gnu gcc/af90  ----------------------
CMACH=PC_LINUX1
F_COMP=af95
F_OPTS=-O3 -YEXT_SFX=_ -YEXT_NAMES=LCS -s -YCFRL=1 -lU77 -OPT:Olimit=0
C_COMP=gcc
C_OPTS=-O3 -DUNDERSCORE -DLITTLE
LOADER=af95
LOADER_OPTS=-O3 -YEXT_SFX=_ -YEXT_NAMES=LCS -s -YCFRL=1 -lU77 -OPT:Olimit=0
LIBS=
#-----------------------------------------------------------------

# For IBM,HP,SGI,ALPHA use these:
#ARCHIVE=ar rs
# For SUN,CONVEX, LINUX
ARCHIVE=ar r
