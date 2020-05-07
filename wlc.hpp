#ifndef WLC
#define WLC

#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
// Utilities and system includes
#include <helper_cuda.h>  // helper function CUDA error checking and initialization
#include <helper_functions.h>  // helper for shared functions common to CUDA Samples

#define WLC_RATIO				(2.8)
#define WLC_STRETCH_CONSTANT	(1.0/WLC_RATIO)
#define WLC_FORCE_CONSTANT		(0.25*WLC_RATIO*WLC_RATIO/((WLC_RATIO-1)*(WLC_RATIO-1))-0.25+1.0/WLC_RATIO)
#define WLC_ENERGY_CONSTANT		(-4.0/(WLC_RATIO*WLC_RATIO)+4.0/WLC_RATIO+1.0/(1-WLC_RATIO))


__device__  double spring_force(double length, double con_length);
__device__ double get_spring_constant(double length, double shear_modulus);

__device__ double get_spring_shear_force(double length, double con_length);


#endif