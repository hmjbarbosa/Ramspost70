#Makefile include include.mk


# RAMS root directory.

#RAMS_ROOT=/path/to/BRAMS32

# Versions.

#BRAMS_VERSION=3.3
#RAMS_VERSION=5.04
#REVU_VERSION=2.3.1
UTILS_VERSION=2.0
#HYPACT_VERSION=1.2.0

# Source directories.

#REVU=$(RAMS_ROOT)/src/post/$(REVU_VERSION)/revu
#COMMON=$(RAMS_ROOT)/src/post/$(REVU_VERSION)/common
#POST_INCS=$(RAMS_ROOT)/src/post/$(REVU_VERSION)/include
#POST_MODS=$(RAMS_ROOT)/src/post/$(REVU_VERSION)/common/modules

#UTILS_LIB=$(RAMS_ROOT)/src/utils/$(UTILS_VERSION)/lib
#EFF=$(RAMS_ROOT)/src/utils/$(UTILS_VERSION)/eff
#NCARGD=$(RAMS_ROOT)/src/utils/$(UTILS_VERSION)/ncargd
#UTILS_INCS=$(RAMS_ROOT)/src/utils/$(UTILS_VERSION)/include
#UTILS_MODS=$(RAMS_ROOT)/src/utils/$(UTILS_VERSION)/lib/modules

#MODEL_LIB=$(RAMS_ROOT)/src/brams/$(BRAMS_VERSION)/lib
#MODEL=$(RAMS_ROOT)/src/brams/$(BRAMS_VERSION)/model
#MODEL_MODS=$(RAMS_ROOT)/src/brams/$(BRAMS_VERSION)/modules
#ISAN=$(RAMS_ROOT)/src/brams/$(BRAMS_VERSION)/isan
#ISAN_MODS=$(RAMS_ROOT)/src/brams/$(BRAMS_VERSION)/modules
#MODEL_INCS=$(RAMS_ROOT)/src/brams/$(BRAMS_VERSION)/include
RCIO=./lib_sources
RCIO_INC=../include/

# Needed for Shallow Cumulus Param.
#SHALLOW=$(RAMS_ROOT)/src/rams/$(RAMS_VERSION)/braz_modules/shallow_cum

# Needed for Grell Cumulus Param.
#GRELL=$(RAMS_ROOT)/src/rams/$(RAMS_VERSION)/braz_modules/grell

# Needed for Heterogenous Soil Moisture Init.
#SOIL=$(RAMS_ROOT)/src/rams/$(RAMS_VERSION)/braz_modules/soil_moisture_init

# Needed for SIB
#SIB=$(RAMS_ROOT)/src/rams/$(RAMS_VERSION)/braz_modules/sib
#New SiB for BRAMS31
#SIB_BRAMS31=$(RAMS_ROOT)/src/rams/$(RAMS_VERSION)/braz_modules/sib_brams31

#For optmizations
#OPTIM=$(RAMS_ROOT)/src/rams/$(RAMS_VERSION)/braz_modules/optim

#For Instrumentation
#RELOGIO=$(RAMS_ROOT)/src/brams/$(BRAMS_VERSION)/relogio


#HYPACT=$(RAMS_ROOT)/src/hypact/$(HYPACT_VERSION)/model
#HYPACT_INCS=$(RAMS_ROOT)/src/hypact/$(HYPACT_VERSION)/include
#HYPACT_MODS=$(RAMS_ROOT)/src/hypact/$(HYPACT_VERSION)/modules
