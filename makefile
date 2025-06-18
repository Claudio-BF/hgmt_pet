CC=gcc
CXX=g++
BASEFLAGS= -Wall -lm
CFLAGS= $(BASEFLAGS) -O2
CXXFLAGS= $(CFLAGS) -std=c++17 -I/usr/lib/python3.13/site-packages/torch/include -I/usr/lib/python3.13/site-packages/torch/include/torch/csrc/api/include
CDBGFLAGS= $(BASEFLAGS) -g
DEPS=linear_algebra.h helper_functions.h hgmt_structs.h compton_chain_ordering.h read_write.h model_wrapper.hh
TORCH_LIBS= -L/usr/lib/python3.13/site-packages/torch/lib -Wl,-rpath,/usr/lib/python3.13/site-packages/torch/lib -ltorch -ltorch_cpu -lc10 -lpthread -fopenmp

hgmt_lor_creator: src/hgmt_lor_creator.o src/hgmt_structs.o src/linear_algebra.o src/helper_functions.o src/compton_chain_ordering.o src/read_write.o neural_network/model_wrapper.o
	$(CXX) -o hgmt_lor_creator $^ $(CFLAGS) $(TORCH_LIBS) -fopenmp

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

%.o: %.cc $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

