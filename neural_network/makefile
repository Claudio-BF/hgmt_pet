CC=gcc
CXX=g++
BASEFLAGS= -Wall -lm
CFLAGS= $(BASEFLAGS) -O2
CDBGFLAGS= $(BASEFLAGS) -g

# LibTorch configuration
PYTORCH_PATH = /usr/lib/python3.13/site-packages/torch
CXXFLAGS = -std=c++17 -Wall -O2 -I$(PYTORCH_PATH)/include -I$(PYTORCH_PATH)/include/torch/csrc/api/include -fopenmp
TORCH_LIBS = -L$(PYTORCH_PATH)/lib -Wl,-rpath,$(PYTORCH_PATH)/lib -ltorch -ltorch_cpu -lc10 -lpthread

DEPS=linear_algebra.h hgmt_structs.h helper_functions.h compton_chain_ordering.h read_write.h

# C objects
C_OBJECTS = training_data.o ../src/hgmt_structs.o ../src/linear_algebra.o ../src/helper_functions.o ../src/compton_chain_ordering.o ../src/read_write.o

# C++ objects
CXX_OBJECTS = model_wrapper.o

hgmt_debug: $(C_OBJECTS) $(CXX_OBJECTS)
	$(CXX) -o training_data $^ $(CFLAGS) $(TORCH_LIBS) -fopenmp

# C compilation rule
%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

# C++ compilation rule
%.o: %.cc
	$(CXX) -c -o $@ $< $(CXXFLAGS)

clean:
	rm -f *.o ../src/*.o training_data

.PHONY: clean

