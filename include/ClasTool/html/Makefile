# Makefile for the DOCUMENTATION pages.
#
# Broken out of each of the directories, since it became too cumbersome to 
# make indexes that way.

all: docs

help:
	@echo '#############################################'
	@echo '#   Makefile to Create html documentation   #'
	@echo '#############################################'
	@echo 'Did you setup ROOT with: use root/2.22 (or higer)?'

docs:
	root.exe Make_HTML.C

clean:
	rm -rf *.html *.ps *.pdf *.h src html examples *.root

distclean:
	rm -rf *.html *.ps *.pdf *.h src html examples *.root htmldoc incl inh inhmem lib ROOT.css ROOT.js LibraryDependencies.gif


