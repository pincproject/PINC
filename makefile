##
## @file		makefile
## @author		Sigvald Marholm <sigvaldm@fys.uio.no>
## @copyright	University of Oslo, Norway
## @brief		PINC makefile.
## @date		10.10.15
##

EXEC	= pinc
CC		= mpicc
CADD	= # Additional CFLAGS accessible from CLI
COPT	= -O3 # Optimization
CFLAGS	=	-std=c11 -Wall\
			-Ilib/iniparser/src\
			-lm -lgsl -lblas -lhdf5 $(COPT) $(CADD)

SDIR	= src
ODIR	= src
HDIR	= src
LDIR	= lib
DDIR	= doc


HEAD_	= pinc.h pusher.h test.h #multigrid.h
SRC_	= main.c io.c aux.c population.c grid.c pusher.c test.c #multigrid.c


OBJ_	= $(SRC_:.c=.o)

HEAD	= $(patsubst %,$(HDIR)/%,$(HEAD_))
SRC		= $(patsubst %,$(SDIR)/%,$(SRC_))
OBJ		= $(patsubst %,$(ODIR)/%,$(OBJ_))

LIBOBJ_	= iniparser/libiniparser.a
LIBHEAD_= iniparser/src/iniparser.h

LIBOBJ = $(patsubst %,$(LDIR)/%,$(LIBOBJ_))
LIBHEAD = $(patsubst %,$(LDIR)/%,$(LIBHEAD_))

all: $(EXEC) cleantestdata doc

test: $(EXEC)

$(EXEC): $(OBJ) $(LIBOBJ)
	@echo "Linking PINC"
	@$(CC) $^ -o $@ $(CFLAGS)
	@echo "PINC is built"

$(ODIR)/%.o: $(SDIR)/%.c $(HEAD)
	@echo "Compiling $<"
	@./check.sh $<
	@$(CC) -c $< -o $@ $(CFLAGS)

$(LDIR)/iniparser/libiniparser.a: $(LIBHEAD)
	@echo "Building iniparser"
	@cd $(LDIR)/iniparser && $(MAKE) > /dev/null 2>&1

$(DDIR)/doxygen/doxyfile.inc: $(DDIR)/doxygen/doxyfile.mk
	@echo INPUT	= ../../$(SDIR) ../../$(HDIR)	> $(DDIR)/doxygen/doxyfile.inc
	@echo FILE_PATTERNS	= $(HEAD_) $(SRC_) >> $(DDIR)/doxygen/doxyfile.inc

doc: $(HEAD) $(SRC) $(DDIR)/doxygen/doxyfile.inc
	@echo "Making documentation (run \"make pdf\" to get pdf)"
	@cd $(DDIR)/doxygen && doxygen doxyfile.mk > /dev/null 2>&1

pdf: doc
	@echo "Making PDF"
	@cd $(DDIR)/latex && $(MAKE) > /dev/null 2>&1

cleandoc:
	@echo "Cleaning documentation"
	@rm -f $(DDIR)/doxygen/doxyfile.inc
	@rm -fr $(DDIR)/html $(DDIR)/latex

cleantestdata:
	@echo "Cleaning test data"
	@rm -f test_*.h5 parsedump.txt

clean: cleandoc cleantestdata
	@echo "Cleaning compilation files (run \"make veryclean\" to clean more)"
	@rm -f *~ $(SDIR)/*.o $(SDIR)/*~ gmon.out

veryclean: clean
	@echo "Cleaning executable and iniparser"
	@rm -f $(EXEC)
	@cd $(LDIR)/iniparser && $(MAKE) veryclean > /dev/null 2>&1
