#ifndef LBFS
#define LBFS

#include "common_kernels.hpp"

#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
// Utilities and system includes
#include <helper_cuda.h>  // helper function CUDA error checking and initialization
#include <helper_functions.h>  // helper for shared functions common to CUDA Samples

extern __constant__ double lattice_weight[15];

__device__ void populate_lattice_macros_uniform(double u_lattice[], double v_lattice[],
	double w_lattice[], double rho_lattice[], double3 cell_1, double3 cell_2,
	double3 interface_node, int i, int neighbour,
	double3 grad_u_1, double3 grad_u_2, double3 grad_v_1, double3  grad_v_2, double3  grad_w_1, double3  grad_w_2, double3  grad_rho_1, double3  grad_rho_2,
	double4 owner_soln, double4 neighbour_soln, int n_cells, double3 cell_normal, double dt, int bcs_rho_type[], int bcs_vel_type[], double4 bcs_macros[]);

__device__ void populate_lattice_macros(double u_lattice[], double v_lattice[],
	double w_lattice[], double rho_lattice[], double3 cell_1, double3 cell_2,
	double3 interface_node, int i, int neighbour,
	double3 grad_u_1, double3 grad_u_2, double3 grad_v_1, double3  grad_v_2, double3  grad_w_1, double3  grad_w_2, double3  grad_rho_1, double3  grad_rho_2,
	double4 owner_soln, double4 neighbour_soln, int n_cells, double3 cell_normal, double dt, int bcs_rho_type[], int bcs_vel_type[], double4 bcs_macros[],
	double d1, double d173);

__device__ void populate_feq(double u_lattice[], double v_lattice[],
	double w_lattice[], double rho_lattice[],
	double feq_lattice[], int k, double pre_conditioned_gamma);

__device__ void populate_feq_lattice(double u_lattice[], double v_lattice[],
	double w_lattice[], double rho_lattice[], double * feq_lattice,
	int k, int face);

__device__ void populate_feq_lattice_explicit(double u_lattice[], double v_lattice[],
	double w_lattice[], double rho_lattice[], double feq_lattice[]);

__device__ void calculate_flux_at_interface(double3 u_interface, double dt, double pre_conditioned_gamma,
	double rho_interface, double4 &cell_flux, int i,
	double* feq_lattice, double3 cell_normal, double interface_area, double* local_fneq, double visc, int testcase, double visc_factor);

__global__ void update_unstructured_bcs(int n_bc_cells, int n_neighbours, int n_cells, int* mesh_owner, int* bcs_rho_type, int* bcs_vel_type, double4* input_soln,
	double4* bc_arr, double3 * cell_centroid, double channel_diameter);

__global__ void get_bc_gradients(int n_bc_cells, int n_neighbours, int n_cells, int* mesh_owner, int* bcs_rho_type, int* bcs_vel_type, double4* input_soln,
	double3* face_normal_arr, double3* centroid, double4* bc_arr,
	double3* grad_rho_arr, double3* grad_u_arr, double3* grad_v_arr, double3* grad_w_arr);
__device__ void add_LS_contributions(double3 &RHS_rho, double3 &RHS_u, double3 &RHS_v, double3 &RHS_w, double4* src, int i1, int i, int i_nb, double3* RHS_array);
__global__ void time_integration(int n_cells, int rk, int rk_max_t, double* delta_t_local, double4* soln_t0, double4* soln_t1, double4* soln,
	double* res_rho, double* res_u, double* res_v, double* res_w, double delta_t, int gpu_time_stepping);
__global__ void	calc_face_flux(int n, double4* input_soln, double* cell_volume, double* surface_area, int* mesh_owner, int* mesh_neighbour, double3* centroid, double3* face_centroid, double3* face_normal,
	double* streaming_dt,
	double3* grad_rho_arr, double3* grad_u_arr, double3* grad_v_arr, double3* grad_w_arr, int n_cells, double pre_conditioned_gamma, double* local_fneq, double visc,
	double* res_rho, double* res_u, double* res_v, double* res_w, double4* res_face,
	int bcs_rho_type[], int bcs_vel_type[], double4 bcs_arr[], double PI, int testcase,  double * visc_factor, bool preconditioning);
__global__ void
calc_face_flux_bc(int n, double4* input_soln, double* cell_volume, double* surface_area, int* mesh_owner, int* mesh_neighbour, double3* centroid, double3* face_centroid, double3* face_normal,
	double* streaming_dt,
	double3* grad_rho_arr, double3* grad_u_arr, double3* grad_v_arr, double3* grad_w_arr, int n_cells, double pre_conditioned_gamma, double* local_fneq, double visc,
	double* res_rho, double* res_u, double* res_v, double* res_w, double4* res_face,
	int bcs_rho_type[], int bcs_vel_type[], double4 bcs_arr[], double PI, int start_index, int testcase, bool preconditioning);
__global__ void calc_total_residual(int n_cells, double4 res_tot, double* res_rho, double* res_u, double* res_v, double* res_w);
__global__ void check_error(int n_cells, double4* soln);
__global__ void get_cfl_device(int n, double4* input_soln, double* cell_volume, double* delta_t_local, double3* cfl_areas, double factor,
	double max_velocity, double pre_conditioned_gamma, double visc, int gpu_time_stepping);
__global__ void get_interior_gradients(int n_cells, int* gradient_cells, double4* input_soln, double3* RHS_array,
	double* LHS_xx, double* LHS_xy, double* LHS_xz,
	double* LHS_yx, double* LHS_yy, double* LHS_yz,
	double* LHS_zx, double* LHS_zy, double* LHS_zz,
	double3* grad_rho_arr, double3* grad_u_arr, double3* grad_v_arr, double3* grad_w_arr); 	


__global__ void calc_face_flux_x(int n, double4* input_soln, double* cell_volume, double* surface_area, int* mesh_owner, int* mesh_neighbour, double3* centroid, double3* face_centroid, double3* face_normal,
	double* streaming_dt,
	double3* grad_rho_arr, double3* grad_u_arr, double3* grad_v_arr, double3* grad_w_arr, int n_cells, double pre_conditioned_gamma, double* local_fneq, double visc,
	double* res_rho, double* res_u, double* res_v, double* res_w, double4* res_face,
	int bcs_rho_type[], int bcs_vel_type[], double4 bcs_arr[], double PI ,double * visc_factor, int testcase);

__global__ void
calc_face_flux_y(int n, double4* input_soln, double* cell_volume, double* surface_area, int* mesh_owner, int* mesh_neighbour, double3* centroid, double3* face_centroid, double3* face_normal,
	double* streaming_dt,
	double3* grad_rho_arr, double3* grad_u_arr, double3* grad_v_arr, double3* grad_w_arr, int n_cells, double pre_conditioned_gamma, double* local_fneq, double visc,
	double* res_rho, double* res_u, double* res_v, double* res_w, double4* res_face,
	int bcs_rho_type[], int bcs_vel_type[], double4 bcs_arr[], double PI, int start_index, double * visc_factor, int testcase);

__global__ void
calc_face_flux_z(int n, double4* input_soln, double* cell_volume, double* surface_area, int* mesh_owner, int* mesh_neighbour, double3* centroid, double3* face_centroid, double3* face_normal,
	double* streaming_dt,
	double3* grad_rho_arr, double3* grad_u_arr, double3* grad_v_arr, double3* grad_w_arr, int n_cells, double pre_conditioned_gamma, double* local_fneq, double visc,
	double* res_rho, double* res_u, double* res_v, double* res_w, double4* res_face,
	int bcs_rho_type[], int bcs_vel_type[], double4 bcs_arr[], double PI, int start_index, double * visc_factor, int testcase);


__global__ void add_face_flux_to_cell(int n_cells, double4 * res, double * cell_volume, int* mesh_owner, int* mesh_neighbour, double* res_rho, double* res_u, double* res_v, double* res_w);

//owner is always west and normal is west to east
__device__ void populate_lattice_macros_uniform_x(double u_lattice[], double v_lattice[],
	double w_lattice[], double rho_lattice[], double cell_1, double cell_2,
	double interface_node, int i, int neighbour,
	double3 grad_u_1, double3 grad_u_2, double3 grad_v_1, double3  grad_v_2, double3  grad_w_1, double3  grad_w_2, double3  grad_rho_1, double3  grad_rho_2,
	double4 owner_soln, double4 neighbour_soln, double dt);

__device__ void populate_lattice_macros_uniform_y(double u_lattice[], double v_lattice[],
	double w_lattice[], double rho_lattice[], double cell_1, double cell,
	double interface_node, int i, int neighbour,
	double3 grad_u_1, double3 grad_u_2, double3 grad_v_1, double3  grad_v_2, double3  grad_w_1, double3  grad_w_2, double3  grad_rho_1, double3  grad_rho_2,
	double4 owner_soln, double4 neighbour_soln, double dt);

//owner is always west and normal is west to east
__device__ void populate_lattice_macros_uniform_z(double u_lattice[], double v_lattice[],
	double w_lattice[], double rho_lattice[], double cell_1, double cell,
	double interface_node, int i, int neighbour,
	double3 grad_u_1, double3 grad_u_2, double3 grad_v_1, double3  grad_v_2, double3  grad_w_1, double3  grad_w_2, double3  grad_rho_1, double3  grad_rho_2,
	double4 owner_soln, double4 neighbour_soln, double dt);

__device__ void calculate_flux_at_interface_x(double3 u_interface, double dt, double pre_conditioned_gamma,
	double rho_interface, double4 &cell_flux, int i,
	double* feq_lattice, double3 cell_normal, double interface_area, double* local_fneq, double visc, double visc_factor, int testcase);

__device__ void calculate_flux_at_interface_y(double3 u_interface, double dt, double pre_conditioned_gamma,
	double rho_interface, double4 &cell_flux, int i,
	double* feq_lattice, double3 cell_normal, double interface_area, double* local_fneq, double visc, double visc_factor, int testcase);


__device__ void calculate_flux_at_interface_z(double3 u_interface, double dt, double pre_conditioned_gamma,
	double rho_interface, double4 &cell_flux, int i,
	double* feq_lattice, double3 cell_normal, double interface_area, double* local_fneq, double visc, double visc_factor, int testcase);

#endif