#
#  Copyright 2015 NVIDIA Corporation
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#

CUDA_PATH       ?= /usr/local/cuda-10.1
HOST_COMPILER ?= g++
NVCC          := $(CUDA_PATH)/bin/nvcc 


################################################################################

# Gencode arguments
#SMS ?=   32 35 50 60 61 70 75
SMS ?=   35 75
ifeq ($(SMS),)
$(info >>> WARNING - no SM architectures have been specified - waiving sample <<<)
SAMPLE_ENABLED := 0
endif


ifeq ($(GENCODE_FLAGS),)
# Generate SASS code for each SM architecture listed in $(SMS)
$(foreach sm,$(SMS),$(eval GENCODE_FLAGS +=  -gencode arch=compute_$(sm),code=sm_$(sm)))

# Generate PTX code from the highest SM architecture in $(SMS) to guarantee forward-compatibility
HIGHEST_SM := $(lastword $(sort $(SMS)))
ifneq ($(HIGHEST_SM),)
GENCODE_FLAGS += -gencode arch=compute_$(HIGHEST_SM),code=compute_$(HIGHEST_SM)
endif
endif

ifeq ($(SAMPLE_ENABLED),0)
EXEC ?= @echo "[@]"
endif
# GENCODE_FLAGS += -v --ptxas-options=-v 
 GENCODE_FLAGS += -maxrregcount=128

################################################################################

CC	= g++
src	= $(wildcard *.cpp *.cu)
src1= $(src:.cpp=.o)
obj	= $(src1:.cu=.o)
#  -I/home/brendan/CGNS/CGNS-3.3.1/src /home/brendan/CGNS/CGNS-3.3.1/src/lib/libcgns.a 
# -O3  -O2 -g
CCFLAGS	=  -O2  -std=c++11 -I/home/Dropbox/PhD/Code/3DLBFS/build/ -I/usr/include/ -I/home/brendan/boost/include/ -I/home/brendan/tecio/teciosrc/ -I/usr/local/cuda-10.1/include/ -I/usr/local/cuda-10.1/samples/common/inc
LDFLAGS	= /home/brendan/boost/lib/libboost_system.a /home/brendan/boost/lib/libboost_filesystem.a  /home/brendan/tecio/teciosrc/libtecio.a  -L/usr/bin/ld/ -lpthread -L/usr/local/cuda-10.1/lib64/ -lcuda -lcudart 
#debug symbols -lineinfo -g -G 
NVCCFLAGS = -lineinfo -std=c++11  -dc -I/home/brendan/boost/include/ -I/usr/local/cuda-10.1/samples/common/inc


myprog: $(obj)
#	 $(CC) -o $@ $^ $(LDFLAGS) 
	 $(NVCC) --gpu-architecture=sm_75 -o $@ $^ $(LDFLAGS) 

%.o : %.cpp
	$(CC) $(CCFLAGS) -o $@ -c $<
	
%.o : %.cu
	$(NVCC) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<
	
.PHONY: clean
clean:
	rm -f $(obj) myprog

#LBFS: $(obj)
#	echo $(CC)
#	$(CC) $(CCFLAGS) $(ACCFLAGS)  $(ACCFLAGS) -o $@  $(LDFLAGS) $<

#clean:
#	$(RM) $(BIN)

#CCFLAGS	= -O3 -Kieee -std=c++11 -Mprof=ccff  -Iinclude -I/home/brendan/boost_1_64_0/prefix/include -I/usr/#local/tecplot/360ex_2018r1/include -I/home/brendan/Eigen -I/home/brendan/CGNS/CGNS-3.3.1/src -I/usr/#include/
