##
##      DUST Makefile
##      Author: Wojciech Jacek Miloch
##      University of Oslo, Norway 
##      2007
##
FILE_EXT = 
ADDCFLAGS ?=	# Additional flags can be set running for instance "make ADDCFLAGS=-DMPI"
##
EXEC=   dust$(FILE_EXT)
##
CC = mpicc
##
##
CFLAGS= -DINSU $(CFLAGS_TRI) -Wall $(ADDCFLAGS)
#-DRESTART
#-Wall 
#-DMPI
#-DRESTART
#-USTATDUST
##CFLAGS= -DCOND 
##-DBEAM
##CFLAGS= -DEPRO
## 
## The flags are used to compile different versions of the program, which 
## allow for studies of the following different cases:
## -DINSU the insulating dust particle (perfect insulator)
## -DCOND the conducting dust particle (perfect conductor)
## -DEPRO the biased probe characteristic
## -DBEAM  additional cold ion beam
##
## Note that -DBEAM can be used all other cases
## It requires extended input file input.txt
## In studies of dust particle 
## the input parameters for the probe
## are not used in the program and may be of any value
#######################################################
##
## FLAGS FOR THE INCLUDED MODIFIED PROGRAM FOR THE POLYGON TRIANGULATION
#
# DEBUG: turn on debugging output
#
# STANDALONE: run as a separate program. read data from file.
#	      If this flag is False, then use the interface procedure 
#	      triangulate_polygon() instead.

CFLAGS_TRI=-UDEBUG -DSTANDALONE -UCLOCK\
	   -I$(inclpath) -L/lib/pa1.1 -g 
inclpath = ./triang
##	
LIBS = -lm 
##
PROBEOBJECTS=	main.o \
		input.o \
		generate.o \
		grid.o \
		accel.o \
		gauss.o \
		diagn.o \
		flux.o \
		dustg.o \
		./fmg/nrutil.o \
		./fmg/fmg.o \
		restart.o \
		photons.o \
		shortcuts.o \
		spherical.o \
##
#TRIANGULATEOBJECTS= 	
#			$(inclpath)/monotone.o \
#			$(inclpath)/construct.o \
#			$(inclpath)/misc.o \
#			$(inclpath)/tri.o \
##
$(EXEC):	$(TRIANGULATEOBJECTS) $(PROBEOBJECTS)
		$(CC) $(CFLAGS_TRI) $(CFLAGS_TRI) -o $(EXEC) $(PROBEOBJECTS) $(TRIANGULATEOBJECTS) $(LIBS)
##
main.o:		main.c const.h funct.h
		$(CC) -c $(CFLAGS) $(CFLAGS_TRI) main.c
##
input.o:	input.c const.h funct.h
		$(CC) -c $(CFLAGS) input.c
##
generate.o:	generate.c const.h funct.h
		$(CC) -c $(CFLAGS) generate.c
##
grid.o:		grid.c const.h funct.h
		$(CC) -c $(CFLAGS) grid.c
##
accel.o:	accel.c const.h funct.h
		$(CC) -c $(CFLAGS) accel.c
##
gauss.o:	gauss.c const.h funct.h
		$(CC) -c $(CFLAGS) gauss.c
##
diagn.o:	diagn.c const.h funct.h
		$(CC) -c $(CFLAGS) diagn.c
##
flux.o:		flux.c const.h funct.h
		$(CC) -c $(CFLAGS) flux.c
##
dustg.o:	dustg.c const.h funct.h
		$(CC) -c $(CFLAGS) dustg.c
##
fmg.o:	./fmg/fmg.c const.h funct.h ./fmg/nrutil.h
		$(CC) -c $(CFLAGS) ./fmg/fmg.c
##
nrutil.o:	./fmg/nrutil.c const.h funct.h ./fmg/nrutil.h
		$(CC) -c $(CFLAGS) ./fmg/nrutil.c
##
restart.o:	restart.c const.h funct.h
		$(CC) -c $(CFLAGS) restart.c
##
photons.o:	photons.c const.h funct.h
		$(CC) -c $(CFLAGS) photons.c
##
shortcuts.o:	shortcuts.c const.h funct.h
		$(CC) -c $(CFLAGS) shortcuts.c
spherical.o:	spherical.c const.h funct.h
		$(CC) -c $(CFLAGS) spherical.c
#########TRIANGULATION
#monotone.o:	$(inclpath)/monotone.c $(inclpath)/triangulate.h
#		$(CC) -c $(CFLAGS_TRI) $(CFLAGS) $(inclpath)/monotone.c
#misc.o:		$(inclpath)/misc.c $(inclpath)/triangulate.h
#		$(CC) -c $(CFLAGS_TRI) $(inclpath)/misc.c
#construct.o:	$(inclpath)/construct.c $(inclpath)/triangulate.h
#		$(CC) -c $(CFLAGS_TRI) $(inclpath)/construct.c
#tri.o:		$(inclpath)/tri.c $(inclpath)/triangulate.h const.h
#		$(CC) -c $(CFLAGS_TRI) $(inclpath)/tri.c

clean:
	rm -f *~ \#* $(PROBEOBJECTS) $(EXEC) $(TRIANGULATEOBJECTS)
