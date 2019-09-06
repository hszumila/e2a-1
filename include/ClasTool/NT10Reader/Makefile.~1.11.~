#   Makefile for CLASTOOL code.
#  
# Notice:
#
# At CEBAF you must execute the following before this will compile and link.
#
# SETUP YOUR CLAS ENVIRONMENT (source /group/clas/builds/../packages/cms/jlab.csh )
# SETUP YOUR ROOT ENVIRONMENT (use root       )
#
#
# REWRITE:
#
#	Maurik.Holtrop@unh.edu	May-12-2004
# 
#
# NOTES:
#
# SRC_CLASS: contains the source files (*.cc) of the CLASSES. A dictionary is made.
# SRC_FILES: contains the source files (*.cc) for the FUNCTIONS and PROGRAMS. No dictionary.
# INC_FILES: contains include files. This is only used to detect if we need to remake the Makefile_depends.
#
# The depends file will only be used if you keep the last include: 'include Makefile_depends'
# This scheme will still try to generate one (in Makefile_bottom).
#
#
#
SRC_CLASS = TNT10Reader.cc
INC_FILES = TNT10Reader.h NT10_commons.h

ifndef CERN_ROOT
	all:
	@echo "CERN_ROOT NOT DEFINED. Skipping NT10Reader "
else
#
# Libraries that must be loaded along with this one.
#
CERNLIB = $(CERN_ROOT)/lib 
OTHERLOADLIBS = -L$(localslib) -lClasBanks -lVirtualReader -L$(CERNLIB) -lpacklib -llapack -lblas -L$(ROOTSYS)/lib  -lHbook 

ifeq ($(findstring Linux,$(OS_NAME)),Linux)
  OTHERLOADLIBS += -ldl -lnsl
endif

include ../Makefile_top

# INCLUDES+= -I$(ROOTSYS)/hbook/inc

#
# We need to override the loader to a FORTRAN aware version...
# For GCC >= 4 we need to use gfortran, otherwise use g77.
# This decision must be made by the shell, since make cannot use dynamic variables in ifdef
#
LD=`echo "__GNUC__" | gcc -E -P - | gawk '/[[:alnum:]]/{if($$1>3){print "gfortran"}else{print "g77"}}'`
#
CREATED_SHLIB=libNT10Reader.$(DllSuf)
CREATED_LIB=libNT10Reader.a
#
# Specify the default
#
all: shlib
#
include ../Makefile_bottom
include Makefile_depends
#
endif
