#ifndef SPRING_NETWORK
#define SPRING_NETWORK

#include "common_kernels.hpp"

#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
// Utilities and system includes
#include <helper_cuda.h>  // helper function CUDA error checking and initialization
#include <helper_functions.h>  // helper for shared functions common to CUDA Samples

//__global__ void update_spring_network_forces(int total_nodes, double * force_x, double * force_y, double * force_z, double * x, double * y, double * z,
//	double * x_ref, double * y_ref, double *  z_ref,
//	double stiffness, double radius, double pi, int object_nodes, double * vel_x, double * vel_y, double * vel_z, double delta_t, double depth, double * object_area);

__global__ void update_spring_network_forces(int total_nodes, double * force_x, double * force_y, double * force_z, double * x, double * y, double * z,
	//volume variables
	double volume_modulus, double total_volume, double total_volume_0, int * node_neighbours, int * node_degree_of_freedom, int max_freedom,
	//area variables
	double local_area_modulus, double global_area_modulus, double total_area, double total_area_0, double * tet_area, double *tet_area_0, double * normal_x, double * normal_y, double *normal_z,
	int * node_tri_neighbours,
	//shear variablesx
	double  * contour_length, double * spring_constant,
	//bending variables
	double global_bending_modulus, double MAD, double MAD_0, double membrane_thickness, double local_bending_modulus, double * curvature, double * curvature_0, double * spring_angle, double viscosity,
	double * vel_x, double * vel_y, double * vel_z, double wlc_ratio, int * tri_neighbour_minor, int * node_neighbours_minor, double * spring_pow_constants, double * spring_angle_0
);



__global__ void get_area_energy(int total_object_nodes, double *object_tet_area, double *object_tet_area_0, double *area_energy, double global_area_modulus);

__global__ void  get_bending_energy(int total_object_nodes, double * object_node_curvature, double * object_node_curvature_0, double * bending_energy, double * object_nodal_area_0,
	double local_bending_modulus);

__global__ void get_spring_energy(int num_nodes, int * node_neighbours, int * node_degree_of_freedom, double * x, double * y, double * z, int max_freedom,
	double * contour_length, double wlc_ratio, double * spring_energy, double * spring_constant);

__global__ void update_spring_network_tet_parameters(int total_tet_nodes, double * x, double * y, double * z, 
	int * tet_connectivity, double * object_area, double * tet_volume, double * tet_area, double * normal_x, double * normal_y, double *normal_z);

__device__ double3 get_angle_derivative(double3 normal_1, double3 normal_2, double angle, double3 a, double3 b, double3 c);

__global__ void preprocess_spring_network_spring_parameters(int num_nodes, int * node_neighbours, int * node_degree_of_freedom, double * x, double * y, double * z, int max_freedom,
	double * contour_length, double * spring_constant, double * node_area, double * node_curvature, int *tri_neighbours, double * tet_area, double * normal_x, double * normal_y, double * normal_z,
	double * spring_angle, double shear_modulus, double wlc_ratio, bool initials, double * spring_pow_constants);

__global__ void external_loading(int num_nodes, double force, double * external_loading_factors, double *force_x, double *nodal_area, double total_area);

__device__ double get_spring_constant(double length, double shear_modulus, double wlc_ratio);

__device__ double get_pivkin_spring_constant(double shear_modulus, double wlc_ratio, double length);

__device__ double3 get_angle_derivative_minor(double3 normal_1, double3 normal_3, double3 a, double3 b, double3 c, double3 e);

__global__ void get_spring_contour_length(int num_nodes, int * node_neighbours, int * node_degree_of_freedom, double * x, double * y, double * z, int max_freedom,
	double * contour_length, double wlc_ratio);

__device__ void get_spring_pow_constants(int i, int max_freedom, int n, double length, double shear_modulus, double wlc_ratio, double * spring_constant, double * pow_constant);

__device__ double3 get_area_derivative(double3 a, double3 b, double3 normal);

__device__ double3 get_angle_derivative_minor_atomic(double3 normal_1, double3 normal_2, double angle, double3 a, double3 b, double3 c);
__device__ double3 get_angle_derivative_atomic(double3 normal_1, double3 normal_2, double angle, double3 a, double3 b, double3 c);

__global__ void update_spring_network_tet_parameters_derivatives(int total_tet_nodes, double * x, double * y, double * z, int * tet_connectivity, double * object_area, double * tet_volume, double * tet_area,
	double * normal_x, double * normal_y, double *normal_z, double* d_Volume, double * d_Normal);

__device__ void get_d_normals(double3 a, double3 b, double3 c, double3 normal, double * tet_volume, double *d_Normal, int i, double * d_Volume);
__device__ void get_d_norms(double * d_norm, double *d_normal, double norm, double3 normal);

__global__ void get_spring_forces_atomics(int total_tets, double * force_x, double * force_y, double * force_z, double * x, double * y, double * z,
	//volume variables
	double volume_modulus, double total_volume, double total_volume_0, int * spring_connectivity, int * node_degree_of_freedom, int max_freedom,
	//area variables
	double local_area_modulus, double global_area_modulus, double total_area, double total_area_0, double * tet_area, double *tet_area_0, double * normal_x, double * normal_y, double *normal_z,
	int * node_tri_neighbours,
	//shear variables
	double  * contour_length, double * spring_constant,
	//bending variables
	double global_bending_modulus, double MAD, double MAD_0, double membrane_thickness, double local_bending_modulus, double * curvature, double * curvature_0, double * spring_angle, double viscosity,
	double * vel_x, double * vel_y, double * vel_z, double wlc_ratio, int * tri_neighbours_minor, int * node_neighbours_minor, double * pow_constant, double * spring_angle_0, double shear_modulus,
	bool initials, double spontaneous_angle
);

__global__ void get_contour_lengths(int total_springs, double * x, double * y, double * z, int * spring_connectivity,
	double  * contour_length, double * spring_constant, double wlc_ratio, double shear_modulus, bool pivkin, double reference_length

);


__global__ void get_reference_curvatures(int total_springs, double * x, double * y, double * z,
	int * spring_connectivity, double * spring_area, double *spring_area_0, double * curvature, double * curvature_0, double * spring_angle, double * spring_angle_0, bool reference_sphere, double curvature_ratio);

__global__ void get_tet_forces_atomics(int total_tets, double * force_x, double * force_y, double * force_z, double * x, double * y, double * z,
	//volume variables
	double volume_modulus, double total_volume, double total_volume_0, int * tet_connectivity, int * node_degree_of_freedom, int max_freedom,
	//area variables
	double local_area_modulus, double global_area_modulus, double total_area, double total_area_0, double * tet_area, double *tet_area_0, double * normal_x, double * normal_y, double *normal_z,
	int * node_tri_neighbours,
	//shear variables
	double  * contour_length, double * spring_constant,
	//bending variables
	double global_bending_modulus, double MAD, double MAD_0, double membrane_thickness, double local_bending_modulus, double * curvature, double * curvature_0, double * spring_angle, double viscosity,
	double * vel_x, double * vel_y, double * vel_z, double wlc_ratio, int * tri_neighbours_minor, int * node_neighbours_minor, double * pow_constant, double * spring_angle_0
);
#endif