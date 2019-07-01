##
## @file		makefile
## @brief		PINC makefile.
## @author		Sigvald Marholm <sigvaldm@fys.uio.no>
##

CC		= mpicc
COPT	= -O3
DOPT 	= -O0

CLOCAL = 	-Ilib/iniparser/src\
			-lm -lgsl -lblas -lhdf5 -lfftw3
LLOCAL =	-Ilib/iniparser/src\
			-lm -lgsl -lblas -lhdf5 -lfftw3

-include local.mk

EXEC	= pinc
CADD	= # Additional CFLAGS accessible from CLI
CFLAGS	= -g -std=c11 -Wall $(CLOCAL) $(COPT) $(CADD) # Flags for compiling
DFLAGS 	= -g $(DOPT) -fno-eliminate-unused-debug-symbols -std=c11 -Wall $(CLOCAL) $(CADD) #flags for debugging
LFLAGS	= -g -std=c11 -Wall $(LLOCAL) $(COPT) $(CADD) # Flags for linking

SDIR	= src
ODIR	= src/obj
HDIR	= src
LDIR	= lib
DDIR	= doc
TSDIR	= test
TODIR	= test/obj
THDIR	= test

HEAD_	= core.h io.h aux.h population.h grid.h pusher.h multigrid.h object.h spectral.h units.h
SRC_	= io.c aux.c population.c grid.c pusher.c multigrid.c object.c spectral.c units.c
OBJ_	= $(SRC_:.c=.o)
DOC_	= main.dox

TESTOBJ_= test.o $(SRC_:.c=.test.o)
#TESTOBJ_ = test.o aux.test.o
TESTHEAD_ = test.h

HEAD	= $(patsubst %,$(HDIR)/%,$(HEAD_))
SRC		= $(patsubst %,$(SDIR)/%,$(SRC_))
OBJ		= $(patsubst %,$(ODIR)/%,$(OBJ_))
TESTOBJ	= $(patsubst %,$(TODIR)/%,$(TESTOBJ_))
TESTHEAD= $(patsubst %,$(THDIR)/%,$(TESTHEAD_))



LIBOBJ_	= iniparser/libiniparser.a
LIBHEAD_= iniparser/src/iniparser.h

LIBOBJ = $(patsubst %,$(LDIR)/%,$(LIBOBJ_))
LIBHEAD = $(patsubst %,$(LDIR)/%,$(LIBHEAD_))

all: version $(EXEC) cleantestdata doc

debug: version $(EXEC).debug

local: version $(EXEC).local cleantestdata doc

test: version $(EXEC).test cleantestdata doc
	@echo "Running Unit Tests"
	@echo $(TEST)
	@./$(EXEC) input.ini

$(EXEC).test: $(TODIR)/main.test.o $(OBJ) $(TESTOBJ) $(LIBOBJ)
	@echo "Linking Unit Tests"
	@$(CC) $^ -o $(EXEC) $(LFLAGS)
	@echo "PINC is built"

$(EXEC).local: $(ODIR)/main.local.o $(OBJ) $(LIBOBJ)
	@echo "Linking PINC (using main.local.c)"
	@$(CC) $^ -o $(EXEC) $(LFLAGS)
	@echo "PINC is built"

$(EXEC): $(ODIR)/main.o $(OBJ) $(LIBOBJ)
	@echo "Linking PINC"
	@$(CC) $^ -o $@ $(LFLAGS)
	@echo "PINC is built"

$(EXEC).debug: $(ODIR)/main.o $(OBJ) $(LIBOBJ)
	@echo "Linking PINC"
	@$(CC) $^ -o $@ $(DFLAGS)
	@echo "PINC debug version is built"

$(ODIR)/%.o: $(SDIR)/%.c $(HEAD)
	@echo "Compiling $<"
	@echo $(HEAD) | xargs -n1 ./aux/check.sh
	@mkdir -p $(ODIR)
	@./aux/check.sh $<
	@$(CC) -c $< -o $@ $(CFLAGS)

$(TODIR)/%.o: $(TSDIR)/%.c $(HEAD) $(TESTHEAD)
	@echo "Compiling $<"
	@echo $(TESTHEAD) | xargs -n1 ./aux/check.sh
	@mkdir -p $(TODIR)
	@./aux/check.sh $<
	@$(CC) -c $< -o $@ -Isrc $(CFLAGS)

$(LDIR)/iniparser/libiniparser.a: $(LIBHEAD)
	@echo "Building iniparser"
	@cd $(LDIR)/iniparser && $(MAKE) libiniparser.a > /dev/null 2>&1

.phony: version
version:
	@echo "Embedding git version"
	@echo "#define VERSION \"$(shell git describe --abbrev=4 --dirty --always --tags)\"" > $(SDIR)/version.h

$(DDIR)/doxygen/doxyfile.inc: $(DDIR)/doxygen/doxyfile.mk $(THDIR)/test.h $(TSDIR)/test.c $(DDIR)/doxygen/$(DOC_)
	@echo INPUT	= ../../$(SDIR) ../../$(HDIR) ../../$(TSDIR) ../../$(THDIR) ../../$(DDIR)/doxygen > $(DDIR)/doxygen/doxyfile.inc
	@echo FILE_PATTERNS	= $(HEAD_) $(SRC_) $(DOC_) test.h test.c  >> $(DDIR)/doxygen/doxyfile.inc

doc: $(HEAD) $(SRC) $(DDIR)/doxygen/doxyfile.inc
	@echo "Making documentation (run \"make pdf\" to get pdf)"
	@cd $(DDIR)/doxygen && doxygen doxyfile.mk > /dev/null 2>&1
	@ln -sf doc/html/index.html doc.html

pdf: doc
	@echo "Making PDF"
	cd $(DDIR)/latex && $(MAKE)	# Intentionally verbose to spot LaTeX errors

cleandoc:
	@echo "Cleaning documentation"
	@rm -f $(DDIR)/doxygen/doxyfile.inc
	@rm -fr $(DDIR)/html $(DDIR)/latex
	@rm -f doc.html

cleantestdata:
	@echo "Cleaning test data"
	@mv data/obj.grid.h5 data/temp 2> /dev/null || true
	@rm -f data/*.h5 data/parsedump.txt
	@mv data/temp data/obj.grid.h5 2> /dev/null || true

clean: cleandoc cleantestdata
	@echo "Cleaning compilation files (run \"make veryclean\" to clean more)"
	@rm -f *~ $(TODIR)/*.o $(ODIR)/*.o $(SDIR)/*.o $(SDIR)/*~ gmon.out ut

janclean:
	@echo "Cleaning compilation files (run \"make veryclean\" to clean more)"
	@rm -f *~ $(TODIR)/*.o $(ODIR)/*.o $(SDIR)/*.o $(SDIR)/*~ gmon.out ut

veryclean: clean
	@echo "Cleaning executable and iniparser"
	@rm -f $(EXEC)
	@cd $(LDIR)/iniparser && $(MAKE) veryclean > /dev/null 2>&1
