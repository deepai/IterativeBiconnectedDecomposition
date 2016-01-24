#Non-Cuda include locations
IDIR = -I../include -I./   

#Cuda Cub library include locations
CUB = -I../cub/

#Cuda cub files dependencies
CUB_H = ../cub/*

#Header files dependencies
includes = ../include/*.h ./*.h

#cpp source files
CPP_SRC = *.cpp
MM_SRC  = ../include/mmio.cpp

EX_INCLUDE = ../include


#cuda source files
CU_SRC = *.cu

#Linker flags
LDFLAGS = -L/usr/local/cuda/lib64 -lcudart

CFLAGS  = -Xcompiler -std=c++11

OMP_FLAGS = -Xcompiler -fopenmp -Xcompiler -lgomp -Xcompiler -w

#DFLAGS = -DINFO

_EXTRA_DEPS = mmio.o Host_Timer.o
EXTRA_DEPS = $(patsubst %,$(EX_INCLUDE)/%,$(_EXTRA_DEPS))

_SRC_DEPS = mmio.cpp Host_Timer.cpp
SRC_DEPS = $(patsubst %,$(EX_INCLUDE)/%,$(_SRC_DEPS))

Main: Main.o KernelFilter.o $(EXTRA_DEPS)
			nvcc -arch=sm_35 $(OMP_FLAGS) $(LDFLAGS) -o Main Main.o KernelFilter.o $(EXTRA_DEPS) $(DFLAGS)

Main.o: $(CU_SRC) $(includes) $(CUB_H)
								nvcc -arch=sm_35 $(OMP_FLAGS) $(CU_SRC) -c $(CUB) $(IDIR) -w $(DFLAGS) -o Main.o

KernelFilter.o: $(CPP_SRC) $(includes)
			 nvcc -arch=sm_35 $(OMP_FLAGS) $(CFLAGS) -c $(CPP_SRC) $(IDIR) -w $(DFLAGS) -o KernelFilter.o

$(EX_INCLUDE)/%.o: $(EX_INCLUDE)/%.cpp $(DEPS)
	nvcc -arch=sm_35 $(CFLAGS) $(LIBS) $(OMP_FLAGS) $(DFLAGS) -c -o $@ $< 

clean: 
			rm -f *~ core $(INCDIR)/*~  *.o ../include/*.o
