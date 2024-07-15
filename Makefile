CUFLAGS= -g -D_FORCE_INLINES -Xcompiler -fopenmp -lineinfo -m64 -arch sm_60 --ptxas-options=-v -I ../include -I /usr/include/nvidia/gdk/
# --ptxas-options=-v -arch sm_20 #-use_fast_math -m64 
#-maxrregcount 31 

#CUFLAGS= -g -G --ptxas-options=-v -arch sm_20 -m64 # -maxrregcount 31 

#CUFLAGS= -Xcompiler -gdwarf-2 -Xcompiler -g3 -G --ptxas-options=-v -arch sm_20 -m64 
# -maxrregcount 31 
NVML_LIB := /usr/include/nvidia/gdk/

NVML_LIB_L := $(addprefix -L , $(NVML_LIB))

LDFLAGS := -lnvidia-ml $(NVML_LIB_L)

LIBS = -lcudart -lcuda -lcurand -lcublas -lnvidia-ml $(NVML_LIB_L) -lgomp
	
OBJS = simple_time.o

TARGETS	= FINSYS

all: $(TARGETS)

$(TARGETS): $(OBJS) 
	g++ -g $(CFLAGS) -o $@ $^

*.o:*.cpp
	g++ -g $(CFLAGS) -c $<

ReflectCoeff_2.o:ReflectCoeff_2.cu
tree2vector.o:tree2vector.cu
kd_struct.o:kd_struct.cu
raystrace.o:raystrace.cu
integral_gpu.o:integral_gpu.cu
virtualface_gpu.o:virtualface_gpu.cu
cudaMallocFree.o:cudaMallocFree.cu
DivRayTube.o:DivRayTube.cu
GPUWatchThreadFun.o:GPUWatchThreadFun.cu
calcThreadFunction.o:calcThreadFunction.cu
socketFunctions.o:socketFunctions.cu
recvThreadFunction.o:recvThreadFunction.cu
scalfuc.o:scalfuc.cu
martixMulti.o:martixMulti.cu
main.o:main.cu
	nvcc $(CUFLAGS) -g -c tree2vector.cu kd_struct.cu raystrace.cu integral_gpu.cu virtualface_gpu.cu cudaMallocFree.cu DivRayTube.cu main.cu  GPUWatchThreadFun.cu calcThreadFunction.cu socketFunctions.cu recvThreadFunction.cu scalfuc.cu ReflectCoeff_2.cu martixMulti.cu

$(TARGETS): tree2vector.o kd_struct.o raystrace.o integral_gpu.o virtualface_gpu.o cudaMallocFree.o DivRayTube.o main.o socketFunctions.o recvThreadFunction.o scalfuc.o ReflectCoeff_2.o
	nvcc $(CUFLAGS) -g -o $@ $(OBJS) tree2vector.o kd_struct.o raystrace.o integral_gpu.o virtualface_gpu.o cudaMallocFree.o DivRayTube.o main.o GPUWatchThreadFun.o calcThreadFunction.o recvThreadFunction.o socketFunctions.o scalfuc.o ReflectCoeff_2.o martixMulti.o $(LIBS)

run:
	./FINSYS
	
clean:
	
	rm -fr *.mod *.o FINSYS
