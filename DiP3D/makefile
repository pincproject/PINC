##
## @file		makefile
## @author		Sigvald Marholm <sigvaldm@fys.uio.no>
## @copyright	University of Oslo, Norway
## @brief		PINC makefile.
## @date		10.10.15
##
EXEC	= pinc
CC		= gcc
CFLAGS	= -std=c11 -Wall\
		-Isrc \
		-Ilib/iniparser/src \
		-lm
DEPS	= src/pinc.h
OBJ		= src/main.o src/input.o src/aux.o
LIBOBJ	= lib/iniparser/libiniparser.a
LIBDEPS	= lib/iniparser/src/iniparser.h

src/%.o: src/%.c $(DEPS)
	$(CC) -c $< -o $@ $(CFLAGS)

$(EXEC): $(OBJ) $(LIBOBJ)
	$(CC) $^ -o $@ $(CFLAGS)

lib/iniparser/libiniparser.a: $(LIBDEPS)
	cd lib/iniparser && $(MAKE)

clean:
	rm -f $(EXEC) *~ src/*.o src/*~
