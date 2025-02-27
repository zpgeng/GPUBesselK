# Compiler settings
CC = gcc
NVCC = nvcc
MPICC = mpicc

# Compiler flags
CFLAGS = -Wall -O2 -I./include
CFLAGS += -I$(CUDA_PATH)/include
NVCCFLAGS = -O2 -arch=sm_75
MPICFLAGS = -Wall -O2 -I./include

# CUDA paths (adjust these to match your CUDA installation)
CUDA_PATH = /usr/local/cuda
CUDA_INC_PATH = $(CUDA_PATH)/include
CUDA_LIB_PATH = $(CUDA_PATH)/lib64

# Include and library flags
INCLUDES = -I$(CUDA_PATH)/include -I./include
LFLAGS = -L$(CUDA_PATH)/lib64 -lcudart
NVCC = nvcc
NVCCFLAGS = $(INCLUDES)

# Source files
CU_SOURCES = src/logbesselk.cu
C_SOURCES = eval_bessel.c

# Object files
C_OBJECTS = $(C_SOURCES:.c=.o)
CU_OBJECTS = $(CU_SOURCES:.cu=.o)

# Executable name
EXECUTABLE = refined_besselk

# Default target
all: $(EXECUTABLE)

# Linking the program
$(EXECUTABLE): $(C_OBJECTS) $(CU_OBJECTS)
	$(NVCC) $^ -o $@ $(LIBS)

# Compiling C sources
src/%.o: src/%.c
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

# Compiling CUDA sources
src/%.o: src/%.cu
	$(NVCC) $(NVCCFLAGS) $(INC) -c $< -o $@

# Clean up
clean:
	rm -f $(C_OBJECTS) $(CU_OBJECTS) $(EXECUTABLE)

# Compiler settings
# NVCC = nvcc
# CC = gcc
# CFLAGS = -O3 -Wall
# NVCCFLAGS = -O3 -rdc=true

# # CUDA path (update this to match your CUDA installation)
# CUDA_PATH = /usr/local/cuda

# # Include and library paths
# INCLUDES = -I$(CUDA_PATH)/include -Iinclude
# LFLAGS = -L$(CUDA_PATH)/lib64 -lcudart

# # Source files
# CUDA_SOURCES = src/cuda_zcmg.cu src/logbesselk.cu
# C_SOURCES = cov_gen.c
# HEADERS = include/logbesselk.h

# # Object files
# CUDA_OBJECTS = $(CUDA_SOURCES:src/%.cu=obj/%.o)
# C_OBJECTS = $(C_SOURCES:%.c=obj/%.o)

# # Executable name
# EXECUTABLE = covariance_matrix

# # Create obj directory
# $(shell mkdir -p obj)

# # Default target
# all: $(EXECUTABLE)

# # Linking
# $(EXECUTABLE): $(CUDA_OBJECTS) $(C_OBJECTS)
# 	$(NVCC) $(LFLAGS) $^ -o $@

# # Compile CUDA source files
# obj/%.o: src/%.cu $(HEADERS)
# 	$(NVCC) $(NVCCFLAGS) $(INCLUDES) -c $< -o $@

# # Compile C source files
# obj/%.o: %.c $(HEADERS)
# 	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# # Clean up
# clean:
# 	rm -rf obj $(EXECUTABLE)

# .PHONY: all clean