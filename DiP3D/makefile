##
##	PINC Makefile
##	Author: Sigvald Marholm
##	University of Oslo, Norway 
##	2015
##
EXEC	= pinc
CC	= gcc
CFLAGS	= -std=c99 \
		-Ilib/iniparser/src

$(EXEC): lib/iniparser/libiniparser.a src/main.c
	$(CC) $(CFLAGS)  src/main.c lib/iniparser/libiniparser.a -o $(EXEC)

clean:
	rm -f $(EXEC) *.o *~
