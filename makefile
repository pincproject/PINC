##
## @file		makefile
## @author		Sigvald Marholm <sigvaldm@fys.uio.no>
## @copyright	University of Oslo, Norway
## @brief		PINC makefile.
## @date		10.10.15
##
EXEC	= pinc
CC		= mpicc
CFLAGS	= -std=c11 -Wall -O3 \
		-Ilib/iniparser/src \
		-lm -lgsl -lblas -static

SDIR	= src
ODIR	= src
HDIR	= src
LDIR	= lib
DDIR	= doc

HEAD_	= pinc.h
SRC_	= main.c io.c aux.c population.c grid.c
OBJ_	= $(SRC_:.c=.o)

HEAD	= $(patsubst %,$(HDIR)/%,$(HEAD_))
SRC		= $(patsubst %,$(SDIR)/%,$(SRC_))
OBJ		= $(patsubst %,$(ODIR)/%,$(OBJ_))

LIBOBJ_	= iniparser/libiniparser.a
LIBHEAD_= iniparser/src/iniparser.h

LIBOBJ = $(patsubst %,$(LDIR)/%,$(LIBOBJ_))
LIBHEAD = $(patsubst %,$(LDIR)/%,$(LIBHEAD_))

all: $(EXEC) doc

$(EXEC): $(OBJ) $(LIBOBJ)
	$(CC) $^ -o $@ $(CFLAGS)
	@echo "PINC is built"

$(ODIR)/%.o: $(SDIR)/%.c $(HEAD)
	./check.sh $<
	$(CC) -c $< -o $@ $(CFLAGS)

$(LDIR)/iniparser/libiniparser.a: $(LIBHEAD)
	cd $(LDIR)/iniparser && $(MAKE)

$(DDIR)/doxyfile.inc: $(DDIR)/doxygen/doxyfile.mk
	@echo INPUT	= ../../$(SDIR) ../../$(HDIR)	> $(DDIR)/doxygen/doxyfile.inc
	@echo FILE_PATTERNS	= $(HEAD_) $(SRC_) >> $(DDIR)/doxygen/doxyfile.inc

doc: $(HEAD) $(SRC) $(DDIR)/doxyfile.inc
	@echo "Making documentation"
	@cd $(DDIR)/doxygen && doxygen doxyfile.mk > /dev/null 2>&1
	@cd $(DDIR)/latex && $(MAKE) > /dev/null 2>&1

cleandoc:
	rm -f $(DDIR)/doxygen/doxyfile.inc
	rm -fr $(DDIR)/html $(DDIR)/latex

clean: cleandoc
	rm -f $(EXEC) *~ $(SDIR)/*.o $(SDIR)/*~

veryclean: clean
	cd $(LDIR)/iniparser && $(MAKE) veryclean
