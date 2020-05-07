#include "wlc.hpp"
#include <cmath>

__device__  double spring_force(double length, double con_length) {

	const double stretchRatio = length / con_length;

	
	
	return ((0.25*pow(1.0 - stretchRatio, -2.0) - 0.25 + stretchRatio - WLC_FORCE_CONSTANT));
}


__device__ double get_spring_constant(double length, double shear_modulus) {

	//equation (by Mingzhu) relating shear modulus to spring constant.
	const double cons_WLC = (0.5*pow(1.0 - WLC_STRETCH_CONSTANT, -3.0) + 1.0);
	
	//equation (by Dao&Fedosov) relating shear modulus to spring constant. //confirmed to be wrong.
	//const REAL cons_WLC = 0.75/((1-r_WLC)*(1-r_WLC))-0.75+4*r_WLC+0.5*r_WLC/pow(1-r_WLC,3.0);
	//const REAL cons_WLC = 0.75/((1-WLC_RATIO)*(1-WLC_RATIO))-0.75+4*WLC_RATIO+0.5*WLC_RATIO/pow(1-WLC_RATIO,3.0);
	//std::cout << "spring_constant: " << ( (4.0*length*SHEARMODULUS)/(SQRT3*cons_WLC) ) << std::endl;
	//getchar();


	return ((4.0*WLC_RATIO*length*shear_modulus) / (sqrt(3.0)*cons_WLC));
}


__device__ double get_spring_shear_force(double length, double con_length) {

	const double stretchRatio = length / con_length;

	printf("WLC_FORCE_CONSTANT %d \n", WLC_FORCE_CONSTANT);
	printf("WLC_spring_force %d \n", (0.25*pow(1.0 - stretchRatio, -2.0) - 0.25 + stretchRatio - WLC_FORCE_CONSTANT));

	return ((0.25*pow(1.0 - stretchRatio, -2.0) - 0.25 + stretchRatio - WLC_FORCE_CONSTANT));
}