CC=gcc
BASEFLAGS= -Wall -lm
CFLAGS= $(BASEFLAGS) -O2
CDBGFLAGS= $(BASEFLAGS) -g
DEPS=vector_ops.h helper_functions.h hgmt_structs.h llist.h compton_chain_ordering.h

hgmt_debug: hgmt_debug.o ../src/vector_ops.o ../src/helper_functions.o ../src/llist.o
	$(CC) -o hgmt_debug $^  $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
