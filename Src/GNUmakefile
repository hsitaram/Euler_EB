# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME = /Users/hsitaram/Documents/ExaCT/BoxLib

DEBUG     = FALSE
USE_MPI   = TRUE
USE_OMP   = FALSE
PROFILE   = FALSE
COMP      = gnu
FCOMP     = gnu
DIM       = 3
PRECISION = DOUBLE
EBASE     = euler

#DEFINES += -DEB_SPHERE
include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

include ./Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules
