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
LIBUTILS=$(UTILS_LIB)/libutils-ia64-$(UTILS_VERSION)-ramspost.a

# Machine-dependent options.
#----------------------- ITANIUM-2 -------------------------------
CMACH=PC_LINUX1

F_COMP=efc
F_OPTS=-O3

C_COMP=ecc
C_OPTS=

LOADER=efc
LOADER_OPTS= -Vaxlib

LIBS=
#-----------------------------------------------------------------

# For IBM,HP,SGI,ALPHA use these:
#ARCHIVE=ar rs
# For SUN,CONVEX, LINUX
ARCHIVE=ar r
