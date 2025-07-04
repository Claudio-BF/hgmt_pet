CC=gcc
CXX=g++

# Common flags
COMMON_FLAGS = -Wall -O3
BASEFLAGS = $(COMMON_FLAGS) -lm
CFLAGS = $(BASEFLAGS)
CDBGFLAGS = $(BASEFLAGS) -g

# LibTorch configuration
PYTORCH_INCLUDES = -I$(PYTORCH_PATH)/include -I$(PYTORCH_PATH)/include/torch/csrc/api/include
CXXFLAGS = $(COMMON_FLAGS) -std=c++17 $(PYTORCH_INCLUDES) -fopenmp
TORCH_LIBS = -L$(PYTORCH_PATH)/lib -Wl,-rpath,$(PYTORCH_PATH)/lib -ltorch -ltorch_cpu -lc10 -lpthread

DEPS=linear_algebra.h hgmt_structs.h helper_functions.h compton_chain_ordering.h read_write.h

# Objects
C_OBJECTS = training_data.o ../src/hgmt_structs.o ../src/linear_algebra.o ../src/helper_functions.o ../src/compton_chain_ordering.o ../src/read_write.o
CXX_OBJECTS = model_wrapper.o

training_data: $(C_OBJECTS) $(CXX_OBJECTS)
	$(CXX) -o $@ $^ $(TORCH_LIBS) -fopenmp -lm

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

%.o: %.cc
	$(CXX) -c -o $@ $< $(CXXFLAGS)
