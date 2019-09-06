#//
#//
#//
#
#
#

ifeq "$(MYSQL_LIBS)" ""
  all:
	@echo PLEASE SETUP MYSQL PATHS: 
	@echo use : setenv MYSQL_LIBS \'lib directory name\'
	@echo and   setenv MYSQL_INCLUDE \'include directory name\'
	@echo ALSO read the readme.
else


ifndef TOP_DIR
  TOP_DIR = $(HOME)
endif
#
#
localobj = $(TOP_DIR)/obj/$(OS_NAME)
localbin = $(TOP_DIR)/bin/$(OS_NAME)
prelocallib = $(TOP_DIR)/lib
locallib = $(TOP_DIR)/lib/$(OS_NAME)

prelocalslib = $(TOP_DIR)/slib
localslib = $(TOP_DIR)/slib/$(OS_NAME)
#
#
CXX      = g++
CPPFLAGS = $(shell root-config --cflags) -I$(MYSQL_INCLUDE)
CXXFLAGS = -g -Wall -fPIC
LDFLAGS  = -g
SOFLAGS  = -shared -Wl,-soname,
INCLUDES = -I../TMySQL
#
#

SOURCE  = TMapTools.cxx
DICT    = TMapToolsDict.cxx
DICTH   = TMapToolsDict.h
DICTO   = $(localobj)/TMapToolsDict.o
OBJECTS = $(localobj)/TMapTools.o $(localobj)/TMapToolsDict.o  $(localobj)/TQuery.o $(localobj)/TQueryDict.o $(localobj)/TItemValue.o $(localobj)/TItemValueDict.o $(localobj)/TItem.o $(localobj)/TItemDict.o $(localobj)/TIntBuffer.o $(localobj)/TIntBufferDict.o $(localobj)/TFloatBuffer.o $(localobj)/TFloatBufferDict.o $(localobj)/TStringBuffer.o  $(localobj)/TStringBufferDict.o


LINKDEF = TMapToolsLinkDef.h

DICTQ     = TQueryDict.cxx
DICTHQ    = TQueryDict.h
DICTOQ    = $(localobj)/TQueryDict.o
OBJECTSQ  = $(localobj)/TQuery.o
LINKDEFQ  = TQueryLinkDef.h
LIB_NAMEQ = libQuery.so

CLASSINCLUDESQ = TQuery.h
CLASSINCLUDES = TMapTools.h
DICTINCLUDES  = -I../TMySQL

LIB_NAME = libMapTools.so

MYSQL = $(MYSQL_LIBS)/libmysqlclient.a

all: $(OBJECTS)
	@test -d $(localslib) || mkdir -p $(localslib)
	$(CXX) $(SOFLAGS) $(LDFLAGS) $(OBJECTS) $(MYSQL) -o $(localslib)/$(LIB_NAME) 

utils: Get_Map_Int Get_Map_Float Scan_Maps Scan_System
	@echo "################################################"
	@echo "#                                              #"
	@echo "#    IMPORTANT: This binaries will apear in    #"
	@echo "#               directory                      #"
	@echo "#  TOP_DIR/bin/OS_NAME                         #"
	@echo "#                                              #"
	@echo "################################################"

Get_Map_Int: $(localobj)/Get_Map_Int.o 
	@test -d $(localbin) || mkdir -p $(localbin)
	$(CXX) $(LDFLAGS) -o $(localbin)/$@ $< $(MYSQL)

Get_Map_Float: $(localobj)/Get_Map_Float.o 
	@test -d $(localbin) || mkdir -p $(localbin)
	$(CXX) $(LDFLAGS) -o $(localbin)/$@ $< $(MYSQL)

Scan_Maps: $(localobj)/Scan_Maps.o 
	@test -d $(localbin) || mkdir -p $(localbin)
	$(CXX) $(LDFLAGS) -o $(localbin)/$@ $< $(MYSQL)

Scan_System: $(localobj)/Scan_System.o 
	@test -d $(localbin) || mkdir -p $(localbin)
	$(CXX) $(LDFLAGS) -o $(localbin)/$@ $< $(MYSQL)

query: $(DICTOQ)  $(OBJECTSQ)
	@test -d $(localslib) || mkdir -p $(localslib)
	$(CXX) $(SOFLAGS) $(LDFLAGS) $(OBJECTSQ) $(DICTOQ) $(MYSQL) -o $(localslib)/$(LIB_NAMEQ) 
	
		

#%.o: %.cxx
#	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -D"__MAKEDUMMLIB__" -c $<

T%Dict.cxx: T%.h
	$(ROOTSYS)/bin/rootcint -f $@ -c $(INCLUDES) $< $(<:.h=LinkDef.h) 


$(localobj)/%.o: %.cxx
	@test -d $(localobj) || mkdir -p $(localobj)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(localobj)/TMapTools*.o core $(localslib)/$(LIB_NAME) TMapToolsDict*.* $(localobj)/TQuery*.o TQueryDict.*  $(localobj)/TMapStorage*.o TMapStorageDict.* $(localobj)/TItemValue*.o TItemValueDict.*

endif