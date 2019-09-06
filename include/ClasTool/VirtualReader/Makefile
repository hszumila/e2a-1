#   Makefile for CLASTOOL code.
#  
# Notice:
#
# At CEBAF you must execute the following before this will compile and link.
#
# SETUP YOUR CLAS ENVIRONMENT (setclas LATEST    )
# SETUP YOUR ROOT ENVIRONMENT (use root/2.22     )
#
#
# Modified:
#      Maurik Holtrop 23/2/2000 :: Restored the help item.
#      Maurik Holtrop 11/4/2000 :: Create ../html , changed CFLAGS for LinuxRH6
#      Maurik Holtrop 28/4/2000 :: Better docs scheme implemented.
#      Maurik Holtrop 30/4/2000 :: Run gcc before rootcint, for better error reporting.

SRC_CLASS = TVirtualReader.cc

include ../Makefile_top

CREATED_LIB = libVirtualReader.a
CREATED_SHLIB = libVirtualReader.$(DllSuf)

OTHERLOADLIBS = -L$(localslib) -lClasBanks
#
#
all: shlib
# docs
include ../Makefile_bottom
#------------- End Of Normal Makefile ------------------------------
#
# Make depend Auto dependency generation.
#
# DO NOT DELETE THIS LINE -- make depend depends on it.
