LIBS = -lm

SRCS=$(wildcard ./src/*.c)

sim: $(SRCS)
	gcc -std=c99 -g -o polar_scl_test main.c $(SRCS) $(LIBS)
	./polar_scl_test

compile: $(SRCS)
	gcc -std=c99 -g -o polar_scl_test main.c $(SRCS) $(LIBS)

compile_ultra: $(SRCS)
	gcc -std=c99 -g -O3 -o polar_scl_test main.c $(SRCS) $(LIBS)
