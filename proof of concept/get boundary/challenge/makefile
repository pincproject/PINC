COPT = -O3

default: main.c
	gcc -o main main.c halo.h halo.o -std=c11 -lm $(COPT)

halo:
	gcc -c halo.c -std=c11 $(COPT)

joint:
	gcc -o main main.c halo.h halo.c -std=c11 -lm $(COPT)

clean:
	rm main
