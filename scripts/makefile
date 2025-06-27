CC=gcc
BASEFLAGS= -Wall -lm
CFLAGS= $(BASEFLAGS) -O2
CDBGFLAGS= $(BASEFLAGS) -g
DEPS=linear_algebra.h helper_functions.h hgmt_structs.h compton_chain_ordering.h read_write.h

hgmt_debug: hgmt_debug.o ../src/linear_algebra.o ../src/hgmt_structs.o ../src/helper_functions.o ../src/read_write.o
	$(CC) -o hgmt_debug $^  $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
