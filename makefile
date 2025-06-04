CC=gcc
BASEFLAGS= -Wall -lm
CFLAGS= $(BASEFLAGS) -O2
CDBGFLAGS= $(BASEFLAGS) -g
DEPS=linear_algebra.h helper_functions.h hgmt_structs.h compton_chain_ordering.h

hgmt_lor_creator: src/hgmt_lor_creator.o src/linear_algebra.o src/helper_functions.o src/compton_chain_ordering.o
	$(CC) -o hgmt_lor_creator $^  $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
