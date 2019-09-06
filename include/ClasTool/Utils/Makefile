#
# The makefile scheme here can really only handle one executable
# per makefile elegantly. Since Utils implies more than one executable,
# we just split the makefile into sub makefiles and call the appropriate one,
# just like in a master makefile.
#

PROGRAMS=WriteRootDst WriteMiniDst Scan_Epics

all: exe

shlib:
	@echo "No library defined."

exe: $(PROGRAMS)

$(PROGRAMS):
	make -f Makefile_$@ exe

dep:
	for i in $(PROGRAMS); do \
		make -f Makefile_$$i dep; \
	done


clean:
	for i in $(PROGRAMS); do \
		make -f Makefile_$$i clean; \
	done

distclean:
	for i in $(PROGRAMS); do \
		make -f Makefile_$$i distclean; \
	done


