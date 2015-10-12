##
## @file		makefile
## @author		Sigvald Marholm <sigvaldm@fys.uio.no>
## @copyright	University of Oslo, Norway
## @brief		PINC makefile.
## @date		10.10.15
##
EXEC	= pinc
CC		= gcc
CFLAGS	= -std=c11 -Wall -O3 \
		-Ilib/iniparser/src \
		-lm

SDIR	= src
ODIR	= src
HDIR	= src
LDIR	= lib
DDIR	= doc

HEAD_	= pinc.h
SRC_	= main.c input.c aux.c
OBJ_	= $(SRC_:.c=.o)

HEAD	= $(patsubst %,$(HDIR)/%,$(HEAD_))
SRC		= $(patsubst %,$(SDIR)/%,$(SRC_))
OBJ		= $(patsubst %,$(ODIR)/%,$(OBJ_))

LIBOBJ_	= iniparser/libiniparser.a
LIBHEAD_= iniparser/src/iniparser.h

LIBOBJ = $(patsubst %,$(LDIR)/%,$(LIBOBJ_))
LIBHEAD = $(patsubst %,$(LDIR)/%,$(LIBHEAD_))

$(EXEC): $(OBJ) $(LIBOBJ)
	$(CC) $^ -o $@ $(CFLAGS)

$(ODIR)/%.o: $(SDIR)/%.c $(HEAD)
	$(CC) -c $< -o $@ $(CFLAGS)

$(LDIR)/iniparser/libiniparser.a: $(LIBHEAD)
	cd $(LDIR)/iniparser && $(MAKE)

$(DDIR)/doxyfile.inc: $(DDIR)/doxyfile.mk
	echo INPUT	= ../$(SDIR) ../$(HDIR)	> $(DDIR)/doxyfile.inc
	echo FILE_PATTERNS	= $(HEAD_) $(SRC_) >> $(DDIR)/doxyfile.inc

doc: $(HEAD) $(SRC) $(DDIR)/doxyfile.inc
	cd $(DDIR) && doxygen doxyfile.mk
	cd $(DDIR)/latex && $(MAKE)

cleandoc:
	rm -f $(DDIR)/doxyfile.inc
	rm -fr $(DDIR)/html $(DDIR)/latex

clean:
	rm -f $(EXEC) *~ $(SDIR)/*.o $(SDIR)/*~
	cd $(LDIR)/iniparser && $(MAKE) veryclean
