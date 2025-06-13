CC=gcc
BASEFLAGS= -Wall -lm
CFLAGS= $(BASEFLAGS) -O2
CDBGFLAGS= $(BASEFLAGS) -g
DEPS=linear_algebra.h hgmt_structs.h helper_functions.h  compton_chain_ordering.h read_write.h 

hgmt_debug: training_data.o ../src/hgmt_structs.o ../src/linear_algebra.o ../src/helper_functions.o ../src/compton_chain_ordering.o ../src/read_write.o
	$(CC) -o training_data $^  $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
