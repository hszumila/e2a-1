#   Makefile for ClasTool
#  
# Notice:
#
# At CEBAF you must execute the following before this will compile and link.
#
# SETUP YOUR CLAS ENVIRONMENT (setclas LATEST    )
# SETUP YOUR ROOT ENVIRONMENT (use root/2.22     )
# LINUX: SETUP EGCS 1.1.1     (use egcs/1.1.1    )
#
# Modified: Maurik Holtrop 4/29/2000  ::
#	Add check for MYSQL_INCLUDE, if not skip those libs.
#	Add target shlib
#       Add targer clean
#

include Makefile_top

ALL_DIRS=ClasBanks VirtualReader DSTReader TClasTool FillBank PartSieve
ALL_DIRS+= MapUtils Filter
SEARCH_INCLUDES=$(ALL_DIRS)
# These require a link to the CLAS programs ($CLAS_PACK) to be available, or they are not 
# needed. Build them with "make extras"
# These are not a required part of the system.
EXTRAS_DIRS = Gsim NT10Reader

EXE_DIRS = Utils

all: include shlib docs

#  shlib: includes clasbanks virtualreader dstreader nt10reader maputils clastool utils

ALL_INCLUDES_LOC=$(foreach dir,$(SEARCH_INCLUDES),$(filter-out %LinkDef.h ,$(wildcard $(dir)/*.h)))
ALL_INCLUDES=$(addprefix include/,$(notdir $(ALL_INCLUDES_LOC)))

help:
	@echo "This is the master makefile for the CLASTOOL package. It calls the individual "
	@echo "Makefiles in each of the subdirectories. "
	@echo "To update the 'include' directory execute 'make include'"
	@echo "To update the dependencies everywhere execute 'make dep'"
	@echo "To cleanup all object and dictionary files execute 'make clean'"
	@echo "To cleanup all libraries as well execute 'make distclean'"
	@echo "To build all the libraries execute 'make'"
	@echo ""
	@echo "I think you are running on a " $(OS_NAME) " system. "
	@echo "so only files for that system will be affected. "
	@echo ""
	@echo "SEARCH_INCLUDES = $(SEARCH_INCLUDES)"
include: $(ALL_INCLUDES_LOC)
	@echo Remaking the include directory....
	@test -d include || mkdir -p include
	@cd include;for i in $(ALL_INCLUDES_LOC); do ln -f -s ../$$i; done

lib: shlib

shlib: 
	for i in $(ALL_DIRS); do \
		make -e -C $$i shlib; \
	done

exe: 
	for i in $(EXE_DIRS); do \
		make -e -C $$i exe; \
	done

extras: 
	for i in $(EXTRAS_DIRS); do \
		make -e -C $$i shlib; \
	done


$(ALL_DIRS):
	@make -e -C $@ all

docs: 
	@make -e -C html docs

dep:
	for i in $(ALL_DIRS); do \
		make -e -C $$i dep; \
	done


clean:
	for i in $(ALL_DIRS) html; do \
		make -e -C $$i clean; \
	done

distclean:
	@rm -f include/*.h
	for i in $(ALL_DIRS) html; do \
		make -e -C $$i distclean; \
	done

.PHONY: include $(ALL_DIRS)


