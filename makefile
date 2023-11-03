VPATH = src/general src/flow_field src/cilia src/cilia/mobility
GEN_FLAGS = -I. -Isrc/general -Isrc/flow_field -Isrc/cilia -Isrc/cilia/mobility -g -w -O3 -lineinfo
CUFCM_ROOT = ../CUFCM/src/

# We only compile the mobility solver that the user has selected.
# This is particularly important when using UAMMD, which takes a long time to compile
# and which the user may not have pulled the source for if uninterested in using it.
MOBILITY_TYPE = $(shell sed -n 's/^ *\#define MOBILITY_TYPE *//p' config.hpp)

ifeq ($(MOBILITY_TYPE), 0) # Stokes drag
MOBILITY_OPTS = -std=c++11
MOBILITY_SOURCE = rpy_mobility_solver.cu stokes_drag_mobility_solver.cu
else ifeq ($(MOBILITY_TYPE), 1) # RPY
MOBILITY_OPTS = -std=c++11
MOBILITY_SOURCE = rpy_mobility_solver.cu
else ifeq ($(MOBILITY_TYPE), 2) # Weakly-coupled-filaments RPY
MOBILITY_OPTS = -std=c++11
MOBILITY_SOURCE = rpy_mobility_solver.cu weakly_coupled_filaments_rpy_mobility_solver.cu
else ifeq ($(MOBILITY_TYPE), 3) # FCM using UAMMD
MOBILITY_OPTS = -std=c++14 -I src/cilia/mobility/UAMMD/src -I src/cilia/mobility/UAMMD/src/third_party -DDOUBLE_PRECISION -lcufft -D_USE_MATH_DEFINES -include ciso646
MOBILITY_SOURCE = 
else ifeq ($(MOBILITY_TYPE), 4) # cuFCM
MOBILITY_OPTS = -arch=sm_75 -std=c++14 -O3 -I../include -lcublas -lcufft -lcblas -lcurand -lcuda -lineinfo
MOBILITY_SOURCE = fcm_mobility_solver.cu $(CUFCM_ROOT)CUFCM_CELLLIST.cu $(CUFCM_ROOT)CUFCM_FCM.cu $(CUFCM_ROOT)CUFCM_DATA.cu $(CUFCM_ROOT)CUFCM_SOLVER.cu $(CUFCM_ROOT)CUFCM_CORRECTION.cu
endif

CILIA_CPP = matrix.cpp quaternion.cpp segment.cpp filament.cpp broyden_solver.cpp rigid_body.cpp swimmer.cpp mobility_solver.cpp
# CILIA_CPP = matrix.cu quaternion.cu util.cu segment.cu filament.cu broyden_solver.cu rigid_body.cu swimmer.cu mobility_solver.cu

CILIA_CUDA = cilia_sim_main.cu seeding.cu cuda_functions.cu globals.cu util.cu $(MOBILITY_SOURCE)

FLOW_FIELD_CPP = flow_field_main.cpp matrix.cpp quaternion.cpp
FLOW_FIELD_CUDA = flow_field_evaluator.cu

cilia_clean:
	-rm cilia
	-rm cilia.exe
	-rm cilia.exp
	-rm cilia.lib
	-rm cilia.pdb

flow_field_clean:
	-rm flow_field
	-rm flow_field.exe
	-rm flow_field.exp
	-rm flow_field.lib
	-rm flow_field.pdb

#
# The following work on the machines at Imperial:
#

# Basic expression works fine on nvidia2.
NVIDIA2_OPTS = -llapack -lblas

cilia_nvidia2: $(CILIA_CPP) $(CILIA_CUDA)
	nvcc $^ $(NVIDIA2_OPTS) $(MOBILITY_OPTS) $(GEN_FLAGS) -o cilia

# Compiles but just outputs zeros if you don't specify the 3.5 architecture of nvidia3's Tesla K40 GPUs.
NVIDIA3_OPTS = -llapack -lblas -arch=sm_35

cilia_nvidia3: $(CILIA_CPP) $(CILIA_CUDA)
	nvcc $^ $(NVIDIA3_OPTS) $(MOBILITY_OPTS) $(GEN_FLAGS) -o cilia

flow_field_nvidia3: $(FLOW_FIELD_CPP) $(FLOW_FIELD_CUDA)
	nvcc $^ $(NVIDIA3_OPTS) $(GEN_FLAGS) -o flow_field

# Basic expression works fine on nvidia4.
NVIDIA4_OPTS = -llapack -lopenblas

cilia_nvidia4: $(CILIA_CPP) $(CILIA_CUDA)
	nvcc $^ $(NVIDIA4_OPTS) $(MOBILITY_OPTS) $(GEN_FLAGS) -o cilia

flow_field_nvidia4: $(FLOW_FIELD_CPP) $(FLOW_FIELD_CUDA)
	nvcc $^ $(NVIDIA4_OPTS) $(GEN_FLAGS) -o flow_field

# On Imperial's HPC cluster, the correct symbolic links don't seem to exist.
# We also have to load the cuda module.
# N.B. Call e.g. "locate liblapack.so" to check for up-to-date versions and their locations if this fails in the future.
HPC_OPTS = -L/usr/lib64 -l:liblapack.so.3.8 -l:libopenblas.so.0

cilia_ic_hpc: $(CILIA_CPP) $(CILIA_CUDA)
	module load cuda/11.4.2 && \
	nvcc $^ $(HPC_OPTS) $(MOBILITY_OPTS) $(GEN_FLAGS) -o cilia

flow_field_ic_hpc: $(FLOW_FIELD_CPP) $(FLOW_FIELD_CUDA)
	module load cuda/11.4.2 && \
	nvcc $^ $(HPC_OPTS) $(GEN_FLAGS) -o flow_field

#
# The following work on my personal PC. Compiling CUDA on Windows requires using the Microsoft
# Visual C++ compiler cl.exe and as such the paths etc. follow the Windows, rather than unix, format.
# blas_win64_MT.dll and lapack_win64_MT.dll have been copied into Cygwin's /bin/
#

PC_OPTS = -ccbin "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Tools\MSVC\14.21.27702\bin\Hostx64\x64" "D:\cygwin64\lib\lib_win64\blas_win64_MT.lib" "D:\cygwin64\lib\lib_win64\lapack_win64_MT.lib"

cilia_pc: $(CILIA_CPP) $(CILIA_CUDA)
	nvcc $^ $(PC_OPTS) $(MOBILITY_OPTS) $(GEN_FLAGS) -o cilia

flow_field_pc: $(FLOW_FIELD_CPP) $(FLOW_FIELD_CUDA)
	nvcc $^ $(PC_OPTS) $(GEN_FLAGS) -o flow_field


# With cuFCM
NVCC_FLAGS=-arch=sm_75 -std=c++17 -O3 -I../include -Xcompiler -fopenmp

LINK=-lcublas -lcufft -llapacke -lcblas -lcurand -lcuda -lineinfo -lopenblas

# CUFCM_FILES = $(CUFCM_ROOT)CUFCM_CELLLIST.cu $(CUFCM_ROOT)CUFCM_FCM.cu $(CUFCM_ROOT)CUFCM_DATA.cu $(CUFCM_ROOT)CUFCM_CORRECTION.cu $(CUFCM_ROOT)CUFCM_SOLVER.cu $(CUFCM_ROOT)CUFCM_RANDOMPACKER.cu

# CUFCM_FILES_SIMPLE = $(CUFCM_ROOT)CUFCM_CELLLIST.cu $(CUFCM_ROOT)CUFCM_FCM.cu $(CUFCM_ROOT)CUFCM_DATA.cu $(CUFCM_ROOT)CUFCM_SOLVER.cu $(CUFCM_ROOT)CUFCM_CORRECTION.cu

# cilia_nvidia4_CUFCM: $(CILIA_CPP) $(CILIA_CUDA)
# 	nvcc $^ -DUSE_DOUBLE_PRECISION $(CFLAGS) $(NVCC_FLAGS) $(NVIDIA4_OPTS) $(LINK) $(GEN_FLAGS) -o bin/cilia

cilia_nvidia4_CUFCM: $(CILIA_CPP) $(CILIA_CUDA)
	nvcc $^ $(CFLAGS) $(NVCC_FLAGS) $(NVIDIA4_OPTS) $(LINK) $(GEN_FLAGS) -o bin/cilia