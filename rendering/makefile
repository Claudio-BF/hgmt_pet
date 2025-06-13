CC=gcc
BASEFLAGS= -Wall -lm
CFLAGS= $(BASEFLAGS) -O2
CDBGFLAGS= $(BASEFLAGS) -g
DEPS=linear_algebra.h hgmt_lor_creator.h helper_functions.h  compton_chain_ordering.h read_write.h imager.h

hgmt_debug: imager.o ../src/hgmt_lor_creator.o ../src/linear_algebra.o ../src/helper_functions.o ../src/compton_chain_ordering.o ../src/read_write.o
	$(CC) -o imager $^  $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
