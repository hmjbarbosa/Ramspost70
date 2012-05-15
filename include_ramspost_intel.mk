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
LIBUTILS=$(UTILS_LIB)/libutils-intel-$(UTILS_VERSION)-ramspost.a

# Machine-dependent options.
#-------------  LINUX Intel v. 8.1 free ifort/icc  ---------------
CMACH=PC_LINUX1
F_COMP=ifort
F_OPTS=-assume byterecl -FR -O1 -static
C_COMP=icc
C_OPTS=-O1 -DUNDERSCORE -DLITTLE
LOADER=ifort 
LOADER_OPTS=-assume byterecl -FR -O1 -static  
LIBS=
#-----------------------------------------------------------------

# For IBM,HP,SGI,ALPHA use these:
#ARCHIVE=ar rs
# For SUN,CONVEX, LINUX
ARCHIVE=ar r
