CC=gcc
BASEFLAGS= -Wall -lm
CFLAGS= $(BASEFLAGS) -O2
CDBGFLAGS= $(BASEFLAGS) -g
DEPS=linear_algebra.h helper_functions.h hgmt_structs.h imager.h

hgmt_debug: imager.o ../src/linear_algebra.o ../src/helper_functions.o
	$(CC) -o imager $^  $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
