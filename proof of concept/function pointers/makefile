default: main.c ../stool/stool.c
	gcc -o main main.c ../stool/stool.c -std=c11 -lm -O0

opt: main.c ../stool/stool.c
	gcc -o main main.c ../stool/stool.c -std=c11 -lm -O3

clean:
	rm main
