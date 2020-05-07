
#include "Solver_gpu.h"
#include <math.h>
#include <cmath>
#include "Solver.h"

#include "vector_var.h"
#include <iostream>
#include "Solution.h"
#include <fstream>
#include "global_variables.h"
#include "residuals.h"
#include <cstdio>
#include <ctime>
#include "artificial_dissipation.h"
#include <boost/math/special_functions/sign.hpp>
#include <limits>
#include "RungeKutta.h"
#include "tecplot_output.h"
#include "gradients.h"
#include <string>
#include <sstream>
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
// Utilities and system includes
#include <helper_cuda.h>  // helper function CUDA error checking and initialization
#include <helper_functions.h>  // helper for shared functions common to CUDA Samples
#include <cuda_profiler_api.h>

#include "lagrangian_object.h"
#include "common_kernels.hpp"
#include "LBFS.hpp"
#include "immersed_boundary_method.hpp"
#include "spring_network.hpp"



using namespace std;

extern __constant__ double lattice_weight[15];


gpu_solver::gpu_solver()
{
	//ctor
}

gpu_solver::~gpu_solver()
{
	//dtor
}



double gpu_solver::min_block(double * input, int n_blocks) {

	double minimum;
	minimum = input[0];
	for (int i = 0; i < n_blocks; i++) {
		minimum = min(input[i], minimum);
	}

	return minimum;
}


double gpu_solver::max_block(double * input, int n_blocks) {

	double maximum;
	maximum = input[0];
	for (int i = 0; i < n_blocks; i++) {
		maximum = max(input[i], maximum);
	}

	return maximum;
}


void gpu_solver::cell_interface_initialiser(double &rho_interface, vector_var &rho_u_interface,
	flux_var &x_flux, flux_var &y_flux) {
	// initialise variables
	 // add in reset function
	rho_interface = 0;

	rho_u_interface.x = 0;
	rho_u_interface.y = 0;
	rho_u_interface.z = 0;

	x_flux.P = 0;
	x_flux.momentum_x = 0;
	x_flux.momentum_y = 0;
	x_flux.momentum_z = 0;



	y_flux.P = 0;
	y_flux.momentum_x = 0;
	y_flux.momentum_y = 0;
	y_flux.momentum_z = 0;

}


double gpu_solver::feq_calc_incomp(double weight, vector_var e_alpha, vector_var u_lattice, double u_magnitude,
	double cs, double rho_lattice, double rho_0, int k) {
	double feq;


	feq = e_alpha.Dot_Product(u_lattice) *3.0;
	feq = feq + (pow(e_alpha.Dot_Product(u_lattice), 2) - pow((u_magnitude* cs), 2))
		*4.5;
	feq = feq * weight *rho_0;
	feq = feq + weight * rho_lattice;

	return feq;

}


double gpu_solver::feq_calc(double weight, vector_var e_alpha, vector_var u_lattice, double u_magnitude,
	double cs, double rho_lattice) {
	double feq;


	feq = 1.0;
	feq = feq
		+ e_alpha.Dot_Product(u_lattice) *3.0;
	feq = feq + (pow(e_alpha.Dot_Product(u_lattice), 2) - pow((u_magnitude* cs), 2))
		*4.5;
	feq = feq * weight *rho_lattice;

	return feq;

}


//get CFL numbers for inviscid and viscous matrices
// see what time stepping results
void gpu_solver::populate_cfl_areas(Solution &cfl_areas, unstructured_mesh &Mesh) {

	double area_x, area_y, area_z;
	int face;

	for (int i = 0; i < Mesh.get_n_cells(); i++) {
		area_x = 0;
		area_y = 0;
		area_z = 0;

		// time step condition as per OpenFoam calcs
		for (int f = 0; f < Mesh.gradient_faces[i].size(); f++) {
			face = Mesh.gradient_faces[i][f];

			// eigen values as per Zhaoli guo(2004) - preconditioning

			//method as per Jiri Blasek: CFD Principles and Application Determination of Max time Step

			// need to calulate correct direction of face vector

			area_x = area_x + fabs(Mesh.get_face_i(face)*Mesh.get_face_area(face));
			area_y = area_y + fabs(Mesh.get_face_j(face)*Mesh.get_face_area(face));
			area_z = area_z + fabs(Mesh.get_face_k(face)*Mesh.get_face_area(face));

		}

		cfl_areas.add_u(i, area_x / 2);
		cfl_areas.add_u(i, area_y / 2);
		cfl_areas.add_u(i, area_z / 2);

	}

	return;
}



void gpu_solver::General_Purpose_Solver_mk_i(unstructured_mesh &Mesh, Solution &soln, Boundary_Conditions &bcs,
	external_forces &source, global_variables &globals, domain_geometry &domain,
	initial_conditions &init_conds, unstructured_bcs &quad_bcs_orig, int mg,
	Solution &residual, int fmg, post_processing &pp,  std::vector<lagrangian_object> &object_vec)
{

	///Declarations
	RungeKutta rk4;

	Solution residual_worker(Mesh.get_total_cells()); // stores residuals


	//Solution wall_shear_stress(Mesh.get_n_wall_cells());
	gradients grads(Mesh.get_total_cells());
	Solution cfl_areas(Mesh.get_total_cells());

	/// Declarations and initialisations

	flux_var RK;


	double4 *temp_soln, *soln_t0, *soln_t1;
	double *force_x, *force_y, *force_z;


	//mesh related GPU variables
	double3 *d_cfl_areas;
	double3 *cell_centroid;
	double3 *face_normal;
	double3 *face_centroid;
	double *cell_volume;
	double *surface_area;
	int* gradient_stencil;
	int* mesh_owner;
	int* mesh_neighbour;
	double *streaming_dt;
	double * plateau_x, *plateau_y, *plateau_z;   /// coordinates for plateau function

	//residual related GPU variables
	double *res_rho, *res_u, *res_v, *res_w;

	///gradient related GPU variables
	double3 *RHS_arr;
	double3 *grad_rho_arr;
	double3 *grad_u_arr;
	double3 *grad_v_arr;
	double3 *grad_w_arr;
	double4 *res_face;
	double *LHS_xx;
	double *LHS_xy;
	double *LHS_xz;
	double *LHS_yx;
	double *LHS_yy;
	double *LHS_yz;
	double *LHS_zx;
	double *LHS_zy;
	double *LHS_zz;

	//bcs related GPU variables
	double4 *bcs_arr;
	int* bcs_rho_type;
	int* bcs_vel_type;


	double4* cell_flux_arr;

	double delta_t = globals.time_marching_step;

	double *d_delta_t_local;
	double *local_fneq;
	double * delta_t_local;
	int *delta_t_frequency;


	/// assign memory
	{

		delta_t_local = new double[Mesh.get_n_cells()];
		if (delta_t_local == NULL) exit(1);
		
		temp_soln = new double4[Mesh.get_total_cells()];
		if (temp_soln == NULL) exit(1);
		soln_t0 = new double4[Mesh.get_total_cells()];
		if (soln_t0 == NULL) exit(1);
		soln_t1 = new double4[Mesh.get_total_cells()];
		if (soln_t1 == NULL) exit(1);
		local_fneq = new double[Mesh.get_total_cells()];
		if (local_fneq == NULL) exit(1);

	}
	//lagrangian object allocations
	double * visc_factor;

	// first get total of object nodes for all cells
	//loop through vector
	int total_object_nodes = 0;
	int total_object_springs = 0;
	int total_object_tets = 0;

	for (int i = 0; i < object_vec.size(); i++) {
		total_object_nodes = total_object_nodes + object_vec[i].num_nodes;
		total_object_springs = total_object_springs + object_vec[i].num_springs;
		total_object_tets = total_object_tets + object_vec[i].num_tets;
	}

	double * object_x_ref, *object_y_ref, *object_z_ref;
	double * object_x, *object_y, *object_z;
	double * object_x0, *object_y0, *object_z0;
	double * object_vel_x, *object_vel_y, *object_vel_z;
	double * object_force_x, *object_force_y, *object_force_z;
	int * object_tet_connectivity;
	double * object_tet_volume, *object_tet_area , *object_tet_area_0;
	double *object_nodal_area, *object_nodal_area_0;
	double * tet_normal_x, *tet_normal_y, *tet_normal_z;

	double * spring_contour_lengths, *spring_constants, *spring_angle, *spring_constants_pow, * spring_angle_0;

	double * object_node_curvature, * object_node_curvature_0, *spring_curvature, * spring_curvature_0;
	
	double * node_add_block;
	double * external_loading_factors;

	double *area_energy, *bending_energy, *shear_energy;
	

	//Spring network variables

	int * object_node_neighbours; //neighbouring nodes to node i ; indexing is i*max_degree_of_freedom 
	int * object_node_neighbours_minor; //opposite minor spring node to node i
	int *object_tri_neighbours; // neighbouring tri neighbours
	int *object_tri_neighbours_minor; // neighbouring triangle to triangle in object_tri_neighbours; needed for curvature bending contribution
	int *object_node_freedom;  // number of degrees of freedom at a node
	int *object_spring_connectivity; // node connectivity for each spring ; index is spring k *4;

	double3 mesh_lengths, mesh_origin;

	mesh_lengths.x = domain.X;
	mesh_lengths.y = domain.Y;
	mesh_lengths.z = domain.Z;

	mesh_origin.x = domain.origin_x;
	mesh_origin.y = domain.origin_y;
	mesh_origin.z = domain.origin_z;

	double local_tolerance;

	double* h_lattice_weight;
	h_lattice_weight = new double[15];
	if (h_lattice_weight == NULL) exit(1);

	double time;
	double output_residual_threshold = 0;
	double visc;
	double angular_freq, wom_cos, force;
	double td; // taylor vortex decay time
	double drag_t1; //drag co-efficients

	std::ofstream error_output, vortex_output, max_u, debug_log, energy_output;
	std::string output_dir, decay_dir, max_u_dir;
	output_dir = globals.output_file + "/error.txt";
	vector_var cell_1, cell_2, interface_node, lattice_node, delta_u, delta_v, delta_w, delta_rho;
	vector_var relative_interface;
	vector_var  vel_lattice, rho_u_interface, u_interface;
	vector_var delta_u1, delta_v1, delta_w1, delta_rho1;
	vector_var cell_normal;
	vector_var flux_e_alpha[9];
	vector_var u, v, w, rho;
	std::vector<vector_var> e_alpha;
	std::vector<int> cell_nodes;

	// vector_var flux_e_alpha;
	residuals convergence_residual;
	flux_var x_flux, y_flux, z_flux;
	flux_var cell_flux;

	flux_var debug[4], debug_flux[4], arti_debug[4];
	flux_var dbug[4];
	flux_var int_debug[4];


	int timesteps;
	int wall = 0;
	
	tecplot_output tecplot;

	///Initialisations

	dt = domain.dt; // timestepping for streaming // non-dim equals 1
	c = 1; // assume lattice spacing is equal to streaming timestep
	cs = c / sqrt(3);
	visc = (globals.tau - 0.5) / 3 * domain.dt;


	local_tolerance = globals.tolerance;
	delta_t = globals.time_marching_step;
	timesteps = ceil(globals.simulation_length);
	output_dir = globals.output_file + "/error.txt";
	decay_dir = globals.output_file + "/vortex_error.txt";
	max_u_dir = globals.output_file + "/max_u.txt";
	// error_output.open("/home/brendan/Dropbox/PhD/Test Cases/Couette Flow/error.txt", ios::out);
	error_output.open(output_dir.c_str(), ios::out);
	output_dir = globals.output_file + "/residual_log.txt";
	debug_log.open(output_dir.c_str(), ios::out);
	vortex_output.open(decay_dir.c_str(), ios::out);
	max_u.open(max_u_dir.c_str(), ios::out);
	output_dir = globals.output_file + "/energy_log.txt";
	energy_output.open(output_dir.c_str(), ios::out);

	time = 0;
	angular_freq = visc * pow(globals.womersley_no, 2) / pow(Mesh.get_Y() / 2, 2);
	force = -init_conds.pressure_gradient;

	time = 0;

	td = 100000000000000000;

	grads.pre_fill_LHS_and_RHS_matrix(bcs, Mesh, domain, soln, globals);

	populate_cfl_areas(cfl_areas, Mesh);


	debug_log << "t,rk,i,res_rho,res_u,res_v,res_w,x,y,z, dt,visc,rho,u,v,ux,uy,uz,vx,vy,vz" << endl;

	/// CUDA checks***********************************//////////////////////
	cudaDeviceProp deviceProp;
	int argc;
	const char *argv = " ";
	/*int devID = findCudaDevice(argc, (const char **)argv);

	if (devID < 0) {
		printf("exiting...\n");
		exit(EXIT_SUCCESS);
	}*/

	int  device_count;
	
	checkCudaErrors(cudaGetDeviceCount(&device_count));

	printf(
		"> Core has  %d devices \n\n", device_count);

	checkCudaErrors(cudaGetDeviceProperties(&deviceProp, globals.gpu_device));

	// Statistics about the GPU device
	printf(
		"> GPU device has %d Multi-Processors, SM %d.%d compute capabilities\n\n",
		deviceProp.multiProcessorCount, deviceProp.major, deviceProp.minor);
	
	post_kernel_checks();
	cudaDeviceSynchronize();

	// assign available cuda device
	//find_available_device(globals.gpu_device);
	checkCudaErrors(cudaSetDevice(globals.gpu_device));

	// num bloacks for different gpu kernels
	int blockSize = 128;
	int numBlocks = (Mesh.get_total_cells() + blockSize - 1) / blockSize;
	int n_Cell_Blocks = (Mesh.get_n_cells() + blockSize - 1) / blockSize;
	int n_bc_Blocks = (Mesh.get_num_bc() + blockSize - 1) / blockSize;
	int n_face_Blocks = (Mesh.get_n_faces() + blockSize - 1) / blockSize;
	int n_neighbours_Blocks = (Mesh.get_n_neighbours() + blockSize - 1) / blockSize;
	int n_face_Blocks_x = (Mesh.get_n_face_x() + blockSize - 1) / blockSize;
	int n_face_Blocks_y = (Mesh.get_n_face_y() + blockSize - 1) / blockSize;
	int n_face_Blocks_z = (Mesh.get_n_face_z() + blockSize - 1) / blockSize;

	int n_node_Blocks = (total_object_nodes + blockSize - 1) / blockSize;
	int n_tet_Blocks = (total_object_tets + blockSize - 1) / blockSize;
	int n_spring_Blocks = (total_object_springs +blockSize -1) / blockSize;

	double delta_x = domain.dt * 2;

	double4 convergence;
	convergence.w = 100000000000;

	double *res_rho_block;
	//res_rho_block = new double[n_Cell_Blocks];
	double *res_u_block;
	//res_u_block = new double[n_Cell_Blocks];
	double *res_v_block;
	//res_v_block = new double[n_Cell_Blocks];
	double *res_w_block;
	//res_w_block = new double[n_Cell_Blocks];

	double * min_delta_t_block;
	double * n_cell_block;
	//min_delta_t_block = new double[n_Cell_Blocks];

	double *tet_add_block, * spring_add_block, *spring_area, *spring_area_0;

	double *force_x_block, *force_y_block, *force_z_block;


	//arrrays for CUDA

	{
		checkCudaErrors(cudaMallocManaged(&res_rho_block, n_Cell_Blocks * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&res_u_block, n_Cell_Blocks * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&res_v_block, n_Cell_Blocks * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&res_w_block, n_Cell_Blocks * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&min_delta_t_block, n_Cell_Blocks * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&n_cell_block, n_Cell_Blocks * sizeof(double)));

		checkCudaErrors(cudaMallocManaged(&d_delta_t_local, Mesh.get_total_cells() * sizeof(double)));

		checkCudaErrors(cudaMallocManaged(&d_cfl_areas, Mesh.get_total_cells() * sizeof(double3)));
		checkCudaErrors(cudaMallocManaged(&temp_soln, Mesh.get_total_cells() * sizeof(double4)));
		checkCudaErrors(cudaMallocManaged(&soln_t0, Mesh.get_total_cells() * sizeof(double4)));
		checkCudaErrors(cudaMallocManaged(&soln_t1, Mesh.get_total_cells() * sizeof(double4)));
		checkCudaErrors(cudaMallocManaged(&cell_volume, Mesh.get_total_cells() * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&gradient_stencil, Mesh.get_n_cells() * sizeof(int) * 6));
		checkCudaErrors(cudaMallocManaged(&mesh_owner, Mesh.get_n_faces() * sizeof(int)));
		checkCudaErrors(cudaMallocManaged(&mesh_neighbour, Mesh.get_n_faces() * sizeof(int)));
		checkCudaErrors(cudaMallocManaged(&cell_centroid, Mesh.get_total_cells() * sizeof(double3)));
		checkCudaErrors(cudaMallocManaged(&face_centroid, Mesh.get_n_faces() * sizeof(double3)));
		checkCudaErrors(cudaMallocManaged(&face_normal, Mesh.get_n_faces() * sizeof(double3)));
		checkCudaErrors(cudaMallocManaged(&surface_area, Mesh.get_n_faces() * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&streaming_dt, Mesh.get_n_faces() * sizeof(double)));

		checkCudaErrors(cudaMallocManaged(&plateau_x, (ceil(domain.X / domain.dt*0.5) + 1) * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&plateau_y, (ceil(domain.Y / domain.dt*0.5) + 1) * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&plateau_z, (ceil(domain.Z / domain.dt*0.5) + 1) * sizeof(double)));


		checkCudaErrors(cudaMallocManaged(&cell_flux_arr, Mesh.get_n_faces() * sizeof(double4)));

		checkCudaErrors(cudaMallocManaged(&force_x, Mesh.get_total_cells() * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&force_y, Mesh.get_total_cells() * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&force_z, Mesh.get_total_cells() * sizeof(double)));

		checkCudaErrors(cudaMallocManaged(&res_rho, Mesh.get_total_cells() * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&res_u, Mesh.get_total_cells() * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&res_v, Mesh.get_total_cells() * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&res_w, Mesh.get_total_cells() * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&res_face, Mesh.get_n_faces() * sizeof(double4)));

		checkCudaErrors(cudaMallocManaged(&local_fneq, Mesh.get_total_cells() * sizeof(double)));

		//arrays for bcs
		checkCudaErrors(cudaMallocManaged(&bcs_arr, Mesh.get_num_bc() * sizeof(double4)));
		checkCudaErrors(cudaMallocManaged(&bcs_rho_type, Mesh.get_num_bc() * sizeof(int)));
		checkCudaErrors(cudaMallocManaged(&bcs_vel_type, Mesh.get_num_bc() * sizeof(int)));


		//arrrays for CUDA Gradient
		checkCudaErrors(cudaMallocManaged(&grad_rho_arr, Mesh.get_total_cells() * sizeof(double3)));
		checkCudaErrors(cudaMallocManaged(&grad_u_arr, Mesh.get_total_cells() * sizeof(double3)));
		checkCudaErrors(cudaMallocManaged(&grad_v_arr, Mesh.get_total_cells() * sizeof(double3)));
		checkCudaErrors(cudaMallocManaged(&grad_w_arr, Mesh.get_total_cells() * sizeof(double3)));
		checkCudaErrors(cudaMallocManaged(&RHS_arr, Mesh.get_n_cells() * sizeof(double3) * 6));
		checkCudaErrors(cudaMallocManaged(&LHS_xx, Mesh.get_n_cells() * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&LHS_xy, Mesh.get_n_cells() * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&LHS_xz, Mesh.get_n_cells() * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&LHS_yx, Mesh.get_n_cells() * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&LHS_yy, Mesh.get_n_cells() * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&LHS_yz, Mesh.get_n_cells() * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&LHS_zx, Mesh.get_n_cells() * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&LHS_zy, Mesh.get_n_cells() * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&LHS_zz, Mesh.get_n_cells() * sizeof(double)));

		checkCudaErrors(cudaMallocManaged(&visc_factor, Mesh.get_n_cells() * sizeof(double)));



		//arrays for lagrangian objects

		checkCudaErrors(cudaMallocManaged(&object_tet_connectivity, total_object_tets * 3 * sizeof(int)));
		checkCudaErrors(cudaMallocManaged(&object_nodal_area, total_object_nodes * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&object_nodal_area_0, total_object_nodes * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&object_node_curvature, total_object_nodes * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&object_node_curvature_0, total_object_nodes * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&spring_curvature, total_object_springs * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&spring_curvature_0, total_object_springs * sizeof(double)));


		checkCudaErrors(cudaMallocManaged(&tet_add_block, n_tet_Blocks * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&node_add_block, n_node_Blocks * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&spring_add_block, n_spring_Blocks * sizeof(double)));

		checkCudaErrors(cudaMallocManaged(&force_x_block, n_node_Blocks * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&force_y_block, n_node_Blocks * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&force_z_block, n_node_Blocks * sizeof(double)));

		checkCudaErrors(cudaMallocManaged(&object_node_neighbours, total_object_nodes * object_vec[0].maxFreedom * sizeof(int)));
		checkCudaErrors(cudaMallocManaged(&object_node_neighbours_minor, total_object_nodes * object_vec[0].maxFreedom * sizeof(int)));
		checkCudaErrors(cudaMallocManaged(&object_tri_neighbours, total_object_nodes * object_vec[0].maxFreedom * sizeof(int)));
		checkCudaErrors(cudaMallocManaged(&object_tri_neighbours_minor, total_object_nodes * object_vec[0].maxFreedom * sizeof(int)));
		checkCudaErrors(cudaMallocManaged(&object_node_freedom, total_object_nodes *sizeof(int)));
		checkCudaErrors(cudaMallocManaged(&object_spring_connectivity, total_object_springs * 4 * sizeof(int)));

		checkCudaErrors(cudaMallocManaged(&spring_contour_lengths, total_object_nodes * object_vec[0].maxFreedom * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&spring_constants, total_object_nodes * object_vec[0].maxFreedom * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&spring_constants_pow, total_object_nodes * object_vec[0].maxFreedom * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&spring_angle, total_object_nodes * object_vec[0].maxFreedom * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&spring_angle_0, total_object_springs* sizeof(double)));

		checkCudaErrors(cudaMallocManaged(&spring_area, total_object_springs * 2*sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&spring_area_0, total_object_springs * 2 * sizeof(double)));

		checkCudaErrors(cudaMallocManaged(&area_energy, total_object_springs * sizeof(double))); // store for reduction
		checkCudaErrors(cudaMallocManaged(&bending_energy, total_object_springs * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&shear_energy, total_object_springs * sizeof(double)));
		

		checkCudaErrors(cudaMallocManaged(&object_tet_volume, total_object_tets *  sizeof(double))); // store for reduction
		checkCudaErrors(cudaMallocManaged(&object_tet_area, total_object_tets * sizeof(double))); // store for reduction
		checkCudaErrors(cudaMallocManaged(&object_tet_area_0, total_object_tets * sizeof(double)));

		checkCudaErrors(cudaMallocManaged(&tet_normal_x, total_object_tets * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&tet_normal_y, total_object_tets * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&tet_normal_z, total_object_tets * sizeof(double)));

		checkCudaErrors(cudaMallocManaged(&object_x_ref, total_object_nodes * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&object_y_ref, total_object_nodes * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&object_z_ref, total_object_nodes * sizeof(double)));

		checkCudaErrors(cudaMallocManaged(&object_x, total_object_nodes * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&object_y, total_object_nodes * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&object_z, total_object_nodes * sizeof(double)));

		checkCudaErrors(cudaMallocManaged(&object_x0, total_object_nodes * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&object_y0, total_object_nodes * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&object_z0, total_object_nodes * sizeof(double)));

		checkCudaErrors(cudaMallocManaged(&object_vel_x, total_object_nodes * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&object_vel_y, total_object_nodes * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&object_vel_z, total_object_nodes * sizeof(double)));

		checkCudaErrors(cudaMallocManaged(&object_force_x, total_object_nodes * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&object_force_y, total_object_nodes * sizeof(double)));
		checkCudaErrors(cudaMallocManaged(&object_force_z, total_object_nodes * sizeof(double)));

		checkCudaErrors(cudaMallocManaged(&external_loading_factors, total_object_nodes * sizeof(double)));
		
	}


	populate_e_alpha(e_alpha, h_lattice_weight, c, globals.PI, 15);
	checkCudaErrors(cudaMemcpyToSymbol(lattice_weight, h_lattice_weight, 15 * sizeof(double)));

	/// Sync before CUDA array used
	cudaDeviceSynchronize();
	populate_cfl_areas(d_cfl_areas, Mesh);


	// transfer class members to arrays for CUDA

	{
		soln_to_double(temp_soln, soln, Mesh.get_total_cells());

		mesh_to_array(cell_volume, Mesh, Mesh.get_total_cells(), "volume",globals.testcase);
		mesh_to_array(gradient_stencil, Mesh, Mesh.get_n_cells(), "gradient_stencil", globals.testcase);
		mesh_to_array(mesh_owner, Mesh, Mesh.get_n_faces(), "mesh_owner", globals.testcase);
		mesh_to_array(mesh_neighbour, Mesh, Mesh.get_n_faces(), "mesh_neighbour", globals.testcase);
		mesh_to_array(surface_area, Mesh, Mesh.get_n_faces(), "surface_area", globals.testcase);
		mesh_to_array(streaming_dt, Mesh, Mesh.get_n_faces(), "streaming_dt", globals.testcase);
		mesh_to_array_double(face_normal, Mesh, Mesh.get_n_faces(), "face_normal");
		mesh_to_array_double(cell_centroid, Mesh, Mesh.get_total_cells(), "cell_centroid");
		mesh_to_array_double(face_centroid, Mesh, Mesh.get_n_faces(), "face_centroid");
		mesh_to_array(plateau_x, Mesh, ceil(domain.X / domain.dt*0.5) + 1, "plateau_x", globals.testcase);
		mesh_to_array(plateau_y, Mesh, ceil(domain.Y / domain.dt*0.5) + 1, "plateau_y", globals.testcase);
		mesh_to_array(plateau_z, Mesh, ceil(domain.Z / domain.dt*0.5) + 1, "plateau_z", globals.testcase);

		gradients_to_array(LHS_xx, grads, Mesh.get_n_cells(), "LHS_xx");
		gradients_to_array(LHS_xy, grads, Mesh.get_n_cells(), "LHS_xy");
		gradients_to_array(LHS_xz, grads, Mesh.get_n_cells(), "LHS_xz");
		gradients_to_array(LHS_yx, grads, Mesh.get_n_cells(), "LHS_yx");
		gradients_to_array(LHS_yy, grads, Mesh.get_n_cells(), "LHS_yy");
		gradients_to_array(LHS_yz, grads, Mesh.get_n_cells(), "LHS_yz");
		gradients_to_array(LHS_zx, grads, Mesh.get_n_cells(), "LHS_zx");
		gradients_to_array(LHS_zy, grads, Mesh.get_n_cells(), "LHS_zy");
		gradients_to_array(LHS_zz, grads, Mesh.get_n_cells(), "LHS_zz");
		gradients_to_array_double(RHS_arr, grads, Mesh.get_n_cells(), "RHS_array");

		bcs_to_array_double(bcs_arr, bcs, Mesh.get_num_bc(), "bcs");
		bcs_to_array(bcs_rho_type, bcs, Mesh.get_num_bc(), "rho_type");
		bcs_to_array(bcs_vel_type, bcs, Mesh.get_num_bc(), "vel_type");
	}


	/// Pre time loop operations on spring network model
	if (globals.testcase != 2) {

	
		if (object_vec[0].type == 2) {
			
			lagrangian_spring_network_to_array(object_vec, object_x_ref, object_y_ref, object_z_ref, object_x, object_y, object_z, object_x0, object_y0, object_z0, object_tet_connectivity, object_node_neighbours,
				object_node_freedom, object_spring_connectivity, object_tri_neighbours,external_loading_factors, object_tri_neighbours_minor, object_node_neighbours_minor);
			post_kernel_checks();
			cudaDeviceSynchronize();

			get_contour_lengths << < n_spring_Blocks, blockSize >> > (total_object_springs, object_x, object_y, object_z, object_spring_connectivity,
				spring_contour_lengths, spring_constants, object_vec[0].wlc_ratio, object_vec[0].shear_modulus, object_vec[0].pivkin_length, object_vec[0].reference_length);
			
			// spring curvatures
			if (object_vec[0].mz_importer) {
				object_vec[0].read_initial_mesh_MZ(object_vec[0].curvature_mesh);
			}
			else {
				object_vec[0].read_initial_mesh(object_vec[0].curvature_mesh);
			}

			get_reference_curvatures << < n_spring_Blocks, blockSize >> > (total_object_springs, object_x, object_y, object_z, object_spring_connectivity, spring_area, spring_area_0,
				spring_curvature, spring_curvature_0, spring_angle, spring_angle_0, object_vec[0].sphere_reference, object_vec[0].curvature_ratio);

			//below function is for non-atomic version, faster on older GPU (TESLA) technology

			/*run_spring_network_functions(total_object_nodes, object_force_x, object_force_y, object_force_z, object_x, object_y, object_z,
					object_vec, object_node_neighbours, object_node_freedom, object_tet_area,
					tet_normal_x, tet_normal_y, tet_normal_z, object_tri_neighbours, spring_contour_lengths, spring_constants,
					object_node_curvature, object_node_curvature_0, spring_angle, object_vel_x, object_vel_y, object_vel_z, tet_add_block,
					n_tet_Blocks, blockSize, node_add_block, n_node_Blocks, total_object_tets, object_tet_connectivity,
					object_nodal_area, object_tet_volume, object_tet_area_0, object_tri_neighbours_minor, object_node_neighbours_minor, true,
					spring_constants_pow, spring_angle_0, object_nodal_area_0);*/
	
			post_kernel_checks();
			cudaDeviceSynchronize();

	

			if (globals.restart_analysis == "true") {
				object_vec[0].solution_read(globals);
			}
			else {  
				if (object_vec[0].mz_importer) {
					object_vec[0].read_initial_mesh_MZ(object_vec[0].mesh_file);
				}
				else {
					object_vec[0].read_initial_mesh(object_vec[0].mesh_file);
				}

			}

			if (globals.testcase == 6 || globals.testcase == 7) {
				object_vec[0].apply_optical_tweezers_force();
			
			}

			//send updated values to GPU memory
			lagrangian_spring_network_to_array_restart(object_vec, object_force_x, object_force_y, object_force_z, object_x, object_y, object_z, object_vel_x, object_vel_y, object_vel_z,external_loading_factors);

		}
		else {
			lagrangian_object_to_array(object_vec, object_x_ref, object_y_ref, object_z_ref, object_x, object_y, object_z, object_x0, object_y0, object_z0, object_tet_connectivity);

		}

		//set up initial reference values for area/volume/curvature etc.
		run_spring_network_functions_atomics(total_object_nodes, object_force_x, object_force_y, object_force_z, object_x, object_y, object_z,
			object_vec, object_node_neighbours, object_node_freedom, object_tet_area,
			tet_normal_x, tet_normal_y, tet_normal_z, object_tri_neighbours, spring_contour_lengths, spring_constants,
			spring_curvature, spring_curvature_0, spring_angle, object_vel_x, object_vel_y, object_vel_z, tet_add_block,
			n_tet_Blocks, blockSize, node_add_block, n_node_Blocks, total_object_tets, object_tet_connectivity,
			object_nodal_area, object_tet_volume, object_tet_area_0, object_tri_neighbours_minor, object_node_neighbours_minor, true, spring_constants_pow,
			spring_angle_0, object_nodal_area_0, object_spring_connectivity, total_object_springs, n_spring_Blocks, spring_area, spring_area_0);


	}
	else {

		fill_double << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), visc_factor, 1.0);
		post_kernel_checks();

	}


	cudaProfilerStart();

	clone_a_to_b << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), temp_soln, soln_t1); // soln_t0 holds macro variable solution at start of time step
	
	// loop in time
	for (int t = 0; t < timesteps; t++) {
		// soln_t0 is the solution at the start of every
		// RK step.(rk = n) Temp_soln holds the values at end of
		// step.(rk = n+1)
		clone_a_to_b << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), soln_t1, soln_t0);// soln_t0 holds macro variable solution at start of time step
		post_kernel_checks();


		//womersley flow peculiarities
		if (globals.testcase == 4) {
			wom_cos = cos(angular_freq * t * delta_t);
			force = -init_conds.pressure_gradient * wom_cos;
		}

		//local timestepping calculation
		// can be removed for uniform grids and replaced with a single calc
		
				
		if (globals.gpu_time_stepping != 2) {

			get_cfl_device << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), temp_soln, cell_volume, d_delta_t_local, d_cfl_areas, globals.time_marching_step,
				globals.max_velocity, globals.pre_conditioned_gamma, globals.visc, globals.gpu_time_stepping);
			post_kernel_checks();
			cudaDeviceSynchronize();
			min<128> << < n_Cell_Blocks, blockSize >> > (d_delta_t_local, min_delta_t_block, Mesh.get_n_cells());
			post_kernel_checks();
			
			delta_t = min_block(min_delta_t_block, n_Cell_Blocks);

		}
		
		// min
		if (globals.gpu_time_stepping == 3) {
			fill_double << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), d_delta_t_local, delta_t);
			//constant
		}else if(globals.gpu_time_stepping == 2) {
			fill_double << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), d_delta_t_local, globals.time_marching_step);
			delta_t = globals.time_marching_step;
		}
		post_kernel_checks();
		cudaDeviceSynchronize();

		if (globals.testcase != 2) {

			if (t == 0) {

				if (object_vec[0].type == 2) {
					tecplot.spring_network_output_gpu(globals.output_file, false, object_vel_x, object_vel_y, object_vel_z, object_x, object_y, object_z, object_force_x, object_force_y,
						object_force_z, object_vec[0].name, 0, object_vec[0].num_nodes, object_vec[0].num_tets, object_tet_connectivity, object_nodal_area, false, object_node_curvature);
				}
				else {
					tecplot.tecplot_output_lagrangian_object_gpu(object_vel_x, object_vel_y, object_vel_z, object_x, object_y, object_z, object_force_x, object_force_y, object_force_z, globals, domain, 0
						, object_vec[0].name, object_vec[0].num_nodes, object_vec[0].depth_nodes, object_vec[0].radial_nodes);
				}

				post_kernel_checks();
				cudaDeviceSynchronize();
			}

			if (globals.testcase != 7) {
				//	need to propogate node position based on final RK4 velocity
				if (globals.mesh_type == 4) {
					interpolate_velocities_on_nodes_non_uniform_grid << < n_node_Blocks, blockSize >> > (total_object_nodes, object_vel_x, object_vel_y, object_vel_z, object_x, object_y, object_z,
						mesh_origin, mesh_lengths, delta_x, temp_soln, Mesh.get_n_cells(), globals.PI, cell_centroid, cell_volume, plateau_x, plateau_y, plateau_z);

				}
				else {
					interpolate_velocities_on_nodes << < n_node_Blocks, blockSize >> > (total_object_nodes, object_vel_x, object_vel_y, object_vel_z, object_x, object_y, object_z,
						mesh_origin, mesh_lengths, delta_x, temp_soln, Mesh.get_n_cells());

					//interpolate_velocities_on_nodes_cos_kernel << < n_node_Blocks, blockSize >> > (total_object_nodes, object_vel_x, object_vel_y, object_vel_z, object_x, object_y, object_z,
					//	mesh_origin, mesh_lengths, delta_x, temp_soln, globals.PI);

				}

				/*		interpolate_velocities_on_nodes_cos_kernel << < n_node_Blocks, blockSize >> > (total_object_nodes, object_vel_x, object_vel_y, object_vel_z, object_x, object_y, object_z,
							mesh_origin, mesh_lengths, delta_x, temp_soln, globals.PI);*/

				update_node_positions << < n_node_Blocks, blockSize >> > (total_object_nodes, object_vel_x, object_vel_y, object_vel_z, object_x, object_y, object_z, delta_t, object_vec[0].num_nodes, 1);
			}
			

		}

		/*fill_zero << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), force_x);*/
		fill_double << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), force_x, init_conds.force_x);
		post_kernel_checks();
		fill_double << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), force_y, init_conds.force_y);
		post_kernel_checks();
		fill_double << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), force_z, init_conds.force_z);
		post_kernel_checks();

		fill_double << < n_node_Blocks, blockSize >> > (total_object_nodes, object_force_x, 0.0);
		post_kernel_checks();
		fill_double << < n_node_Blocks, blockSize >> > (total_object_nodes, object_force_y, 0.0);
		post_kernel_checks();
		fill_double << < n_node_Blocks, blockSize >> > (total_object_nodes, object_force_z, 0.0);
		post_kernel_checks();

		if (globals.testcase != 2) {

			if (object_vec[0].type == 2) {

				post_kernel_checks();
				cudaDeviceSynchronize();


				run_spring_network_functions_atomics(total_object_nodes, object_force_x, object_force_y, object_force_z, object_x, object_y, object_z,
					object_vec, object_node_neighbours, object_node_freedom, object_tet_area,
					tet_normal_x, tet_normal_y, tet_normal_z, object_tri_neighbours, spring_contour_lengths, spring_constants,
					spring_curvature, spring_curvature_0, spring_angle, object_vel_x, object_vel_y, object_vel_z, tet_add_block,
					n_tet_Blocks, blockSize, node_add_block, n_node_Blocks, total_object_tets, object_tet_connectivity,
					object_nodal_area, object_tet_volume, object_tet_area_0, object_tri_neighbours_minor, object_node_neighbours_minor, false,
					spring_constants_pow, spring_angle_0, object_nodal_area_0, object_spring_connectivity, total_object_springs, n_spring_Blocks, spring_area, spring_area_0);

				/*run_spring_network_functions(total_object_nodes, object_force_x, object_force_y, object_force_z, object_x, object_y, object_z,
					object_vec, object_node_neighbours, object_node_freedom, object_tet_area,
					tet_normal_x, tet_normal_y, tet_normal_z, object_tri_neighbours, spring_contour_lengths, spring_constants,
					object_node_curvature, object_node_curvature_0, spring_angle, object_vel_x, object_vel_y, object_vel_z, tet_add_block,
					n_tet_Blocks, blockSize, node_add_block, n_node_Blocks, total_object_tets, object_tet_connectivity,
					object_nodal_area, object_tet_volume, object_tet_area_0, object_tri_neighbours_minor, object_node_neighbours_minor, false,
					spring_constants_pow, spring_angle_0, object_nodal_area_0);*/


				post_kernel_checks();
				cudaDeviceSynchronize();

				//optical tweezers, add external forces
				if (globals.testcase == 6 || globals.testcase == 7) {
					external_loading << < n_node_Blocks, blockSize >> > (total_object_nodes, object_vec[0].boundary_force, external_loading_factors, object_force_x, object_nodal_area, object_vec[0].total_area);

				}

				//really shit code, but I'm only using one RBC so it's acceptable to limit atomics.
				// needs updating for multi RBC, probably fine on Volta GPUs
				if (globals.testcase != 7) {
					if (t % 10000 == 1) {
						update_viscosities(object_x, object_y, object_z, n_cell_block, n_Cell_Blocks, n_node_Blocks, Mesh, visc_factor, cell_centroid,
							total_object_nodes, node_add_block, blockSize, Mesh.get_delta_h(), globals.PI, object_vec[0].internal_viscosity_ratio, total_object_tets, object_tet_connectivity);
					}

				}
			}
			else {
				//for now assume uniform stiffness, radius etc. 
				update_node_forces << < n_node_Blocks, blockSize >> > (total_object_nodes, object_force_x, object_force_y, object_force_z, object_x, object_y, object_z, object_x_ref, object_y_ref, object_z_ref,
					object_vec[0].stiffness, object_vec[0].radius, globals.PI, object_vec[0].num_nodes, object_vel_x, object_vel_y, object_vel_z, delta_t, object_vec[0].depth);
			}


			if (globals.mesh_type == 4) {
				if (globals.testcase != 7) {
					spread_forces_on_non_uniform_grid << < n_node_Blocks, blockSize >> > (total_object_nodes, object_force_x, object_force_y, object_force_z, object_x, object_y, object_z,
						mesh_origin, mesh_lengths, delta_x, force_x, force_y, force_z, Mesh.get_n_cells(), globals.PI, cell_centroid, cell_volume,
						plateau_x, plateau_y, plateau_z);
				}
			}
			else {
				//assume uniform grid for now, need moving least squares stencil in the future
				spread_forces_on_structured_grid << < n_node_Blocks, blockSize >> > (total_object_nodes, object_force_x, object_force_y, object_force_z, object_x, object_y, object_z,
					mesh_origin, mesh_lengths, delta_x, force_x, force_y, force_z, Mesh.get_n_cells());

				/*spread_forces_on_structured_grid_cos_kernel << < n_node_Blocks, blockSize >> > (total_object_nodes, object_force_x, object_force_y, object_force_z, object_x, object_y, object_z,
					mesh_origin, mesh_lengths, delta_x, force_x, force_y, force_z, Mesh.get_n_cells(), globals.PI);*/


			}
			/*spread_forces_on_structured_grid_cos_kernel << < n_node_Blocks, blockSize >> > (total_object_nodes, object_force_x, object_force_y, object_force_z, object_x, object_y, object_z,
				mesh_origin, mesh_lengths, delta_x, force_x, force_y, force_z, Mesh.get_n_cells(), globals.PI);*/

			if (globals.testcase == 7) {
				update_node_positions_tweezers << < n_node_Blocks, blockSize >> > (total_object_nodes, object_vel_x, object_vel_y, object_vel_z, object_x, object_y, object_z,
					delta_t, total_object_nodes, object_vec[0].point_mass, object_force_x, object_force_y, object_force_z);
					 
			}

		}
	


		if (globals.testcase != 7) {
			for (int rk = 0; rk < rk4.timesteps; rk++) {

				drag_t1 = 0.0;

				//update temp_soln boundary conditions
				update_unstructured_bcs << < n_bc_Blocks, blockSize >> > (Mesh.get_num_bc(), Mesh.get_n_neighbours(), Mesh.get_n_cells(), mesh_owner, bcs_rho_type, bcs_vel_type, temp_soln, bcs_arr, cell_centroid, domain.Y);

				//set to zeros
				fill_zero << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), res_rho);
				post_kernel_checks();
				fill_zero << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), res_u);
				post_kernel_checks();
				fill_zero << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), res_v);
				post_kernel_checks();
				fill_zero << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), res_w);
				post_kernel_checks();
				fill_zero << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), res_face);
				post_kernel_checks();

				// time2 = clock();
				get_interior_gradients << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), gradient_stencil, temp_soln,
					RHS_arr, LHS_xx, LHS_xy, LHS_xz, LHS_yx, LHS_yy, LHS_yz, LHS_zx, LHS_zy, LHS_zz,
					grad_rho_arr, grad_u_arr, grad_v_arr, grad_w_arr);
				post_kernel_checks();


				//get boundary condition gradients
				get_bc_gradients << < n_bc_Blocks, blockSize >> > (Mesh.get_num_bc(), Mesh.get_n_neighbours(), Mesh.get_n_cells(), mesh_owner, bcs_rho_type, bcs_vel_type, temp_soln,
					face_normal, cell_centroid, bcs_arr,
					grad_rho_arr, grad_u_arr, grad_v_arr, grad_w_arr);
				post_kernel_checks();

				//time3 = clock();

				//std::cout << "CPU Cycles Gradients:" << double(time3 - time2) << std::endl;
				wall = 0;
				// loop through each cell and exclude the ghost cells
				//using n_cells here rather than total_cells

				
				if (globals.mesh_type == 3) {
					calc_face_flux << < n_neighbours_Blocks, blockSize >> > (Mesh.get_n_neighbours(), temp_soln, cell_volume, surface_area, mesh_owner, mesh_neighbour, cell_centroid, face_centroid, face_normal,
						streaming_dt, grad_rho_arr, grad_u_arr, grad_v_arr, grad_w_arr, Mesh.get_n_cells(), (1 / globals.pre_conditioned_gamma), local_fneq, globals.visc,
						res_rho, res_u, res_v, res_w, res_face,
						bcs_rho_type, bcs_vel_type, bcs_arr, globals.PI, globals.testcase,visc_factor, globals.preconditioning_active);
					post_kernel_checks();
					cudaDeviceSynchronize();
				}
				else {
					calc_face_flux_x << < n_face_Blocks_x, blockSize >> > (Mesh.get_n_face_x(), temp_soln, cell_volume, surface_area, mesh_owner, mesh_neighbour, cell_centroid, face_centroid, face_normal,
						streaming_dt, grad_rho_arr, grad_u_arr, grad_v_arr, grad_w_arr, Mesh.get_n_cells(), (1 / globals.pre_conditioned_gamma), local_fneq, globals.visc,
						res_rho, res_u, res_v, res_w, res_face,
						bcs_rho_type, bcs_vel_type, bcs_arr, globals.PI, visc_factor,globals.testcase);
					calc_face_flux_y << < n_face_Blocks_y, blockSize >> > (Mesh.get_n_face_y(), temp_soln, cell_volume, surface_area, mesh_owner, mesh_neighbour, cell_centroid, face_centroid, face_normal,
						streaming_dt, grad_rho_arr, grad_u_arr, grad_v_arr, grad_w_arr, Mesh.get_n_cells(), (1 / globals.pre_conditioned_gamma), local_fneq, globals.visc,
						res_rho, res_u, res_v, res_w, res_face,
						bcs_rho_type, bcs_vel_type, bcs_arr, globals.PI, Mesh.get_n_face_x(), visc_factor,globals.testcase);
					calc_face_flux_z << < n_face_Blocks_z, blockSize >> > (Mesh.get_n_face_z(), temp_soln, cell_volume, surface_area, mesh_owner, mesh_neighbour, cell_centroid, face_centroid, face_normal,
						streaming_dt, grad_rho_arr, grad_u_arr, grad_v_arr, grad_w_arr, Mesh.get_n_cells(), (1 / globals.pre_conditioned_gamma), local_fneq, globals.visc,
						res_rho, res_u, res_v, res_w, res_face,
						bcs_rho_type, bcs_vel_type, bcs_arr, globals.PI, (Mesh.get_n_face_x() + Mesh.get_n_face_y()), visc_factor,globals.testcase);
					post_kernel_checks();
					cudaDeviceSynchronize();

				}

				calc_face_flux_bc << < n_bc_Blocks, blockSize >> > (Mesh.get_num_bc(), temp_soln, cell_volume, surface_area, mesh_owner, mesh_neighbour, cell_centroid, face_centroid, face_normal,
					streaming_dt, grad_rho_arr, grad_u_arr, grad_v_arr, grad_w_arr, Mesh.get_n_cells(), (1 / globals.pre_conditioned_gamma), local_fneq, globals.visc,
					res_rho, res_u, res_v, res_w, res_face,
					bcs_rho_type, bcs_vel_type, bcs_arr, globals.PI, Mesh.get_n_neighbours(), globals.testcase, globals.preconditioning_active);
				post_kernel_checks();
				cudaDeviceSynchronize();


				add_face_flux_to_cell << < n_face_Blocks, blockSize >> > (Mesh.get_n_faces(), res_face, cell_volume, mesh_owner, mesh_neighbour, res_rho, res_u, res_v, res_w);
				post_kernel_checks();
				cudaDeviceSynchronize();
				
				add << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), res_u, force_x);
				add << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), res_v, force_y);
				add << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), res_w, force_z);
				post_kernel_checks();
				cudaDeviceSynchronize();

				//Update  solutions  //update RK values

				time_integration << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), rk, rk4.timesteps, d_delta_t_local, soln_t0, soln_t1, temp_soln,
					res_rho, res_u, res_v, res_w, delta_t, globals.gpu_time_stepping);

				post_kernel_checks();

			}

		}


		//check for periodic condition
		if (t % 5000 == 1 && object_vec[0].periodic) {
			check_periodic_object(blockSize, object_z, node_add_block, n_node_Blocks, total_object_nodes, domain, n_Cell_Blocks);
		}


		//get square of residuals
		square << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), res_rho);
		square << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), res_u);
		square << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), res_v);
		square << < n_Cell_Blocks, blockSize >> > (Mesh.get_n_cells(), res_w);

		//reduce add residuals
		total<128> << < n_Cell_Blocks, blockSize >> > (res_rho, res_rho_block, Mesh.get_n_cells());
		total<128> << < n_Cell_Blocks, blockSize >> > (res_u, res_u_block, Mesh.get_n_cells());
		total<128> << < n_Cell_Blocks, blockSize >> > (res_v, res_v_block, Mesh.get_n_cells());
		total<128> << < n_Cell_Blocks, blockSize >> > (res_w, res_w_block, Mesh.get_n_cells());
		post_kernel_checks();
		cudaDeviceSynchronize();
		convergence_residual.reset();
		convergence_residual.l2_norm_rms_moukallad(globals, res_rho_block, res_u_block, res_v_block, res_w_block, n_Cell_Blocks, Mesh.get_n_cells());

		time = t * delta_t;

		bool stopped = false;


		if ((object_vel_x[object_vec[0].max_index] < 0.001 || object_vel_x[object_vec[0].min_index] > -0.001) && t > globals.output_step) {
			stopped = true;
		}

		if (mg == 0 && t%globals.output_step == 1) {


			 double volume_coefficient;

			volume_coefficient = -1 * object_vec[0].volume_modulus * (object_vec[0].total_volume / object_vec[0].total_volume_0 - 1.0);

			 double area_coefficient;

			 area_coefficient = -1 * object_vec[0].global_area_modulus * (object_vec[0].total_area / object_vec[0].total_area_0 - 1.0);

			 double MAD_coefficient;

			MAD_coefficient = -1 * object_vec[0].global_bending_modulus * (object_vec[0].total_MAD - object_vec[0].total_MAD_0) / object_vec[0].total_area_0;


			double ADE_energy;

			ADE_energy = 0.5 * object_vec[0].global_bending_modulus * (object_vec[0].total_MAD - object_vec[0].total_MAD_0) * (object_vec[0].total_MAD - object_vec[0].total_MAD_0) /
				 object_vec[0].total_area_0;

			double volume_energy;
			volume_energy = 0.5 * object_vec[0].volume_modulus * (object_vec[0].total_volume - object_vec[0].total_volume_0 ) * (object_vec[0].total_volume - object_vec[0].total_volume_0) /
				 object_vec[0].total_volume_0;

			double global_area_energy;

			global_area_energy = 0.5  * object_vec[0].global_area_modulus *
				(object_vec[0].total_area - object_vec[0].total_area_0)*(object_vec[0].total_area - object_vec[0].total_area_0) / object_vec[0].total_area_0;
						
			double local_area_energy;
			local_area_energy = 0.0;

			get_area_energy <<< n_spring_Blocks, blockSize >> > (total_object_springs, spring_area, spring_area_0, area_energy,  object_vec[0].global_area_modulus) ;

			post_kernel_checks();
			cudaDeviceSynchronize();
			
			total<128> << < n_spring_Blocks, blockSize >> > (area_energy, spring_add_block, total_object_springs);
			post_kernel_checks();
			cudaDeviceSynchronize();
			
			for (int q = 0; q < n_tet_Blocks; q++) {
				local_area_energy += spring_add_block[q];
			}

			double local_bending_energy;
			local_bending_energy = 0.0;

			get_bending_energy << < n_spring_Blocks, blockSize >> > (total_object_springs, spring_curvature, spring_curvature_0, bending_energy, spring_area,
				object_vec[0].local_bending_modulus);
			post_kernel_checks();
			cudaDeviceSynchronize();

			total<128> << < n_spring_Blocks, blockSize >> > (bending_energy, spring_add_block, total_object_springs);
			post_kernel_checks();
			cudaDeviceSynchronize();

			for (int q = 0; q < n_tet_Blocks; q++) {
				local_bending_energy += spring_add_block[q];
			}

			double spring_energy;
			spring_energy = 0.0;


			get_spring_energy << < n_spring_Blocks, blockSize >> >  (total_object_springs, object_spring_connectivity, object_spring_connectivity, object_x, object_y, object_z, object_vec[0].maxFreedom,
				spring_contour_lengths, (1/object_vec[0].wlc_ratio),  shear_energy, spring_constants);

			post_kernel_checks();
			cudaDeviceSynchronize();

			total<128> << < n_spring_Blocks, blockSize >> > (shear_energy, spring_add_block, total_object_springs);
			post_kernel_checks();
			cudaDeviceSynchronize();

			for (int q = 0; q < n_tet_Blocks; q++) {
				spring_energy += spring_add_block[q];
			}

			double total_energy;

			total_energy = volume_energy + global_area_energy + local_area_energy + local_bending_energy + ADE_energy+ spring_energy;

			energy_output <<  t << ", " << volume_energy << ", " <<
				global_area_energy << ", " << local_area_energy << ", " <<
				local_bending_energy << ", " << ADE_energy << ", " <<
				spring_energy << " ,  " << total_energy << endl;

			if (object_vec[0].type == 2) {
				tecplot.spring_network_output_gpu(globals.output_file, false, object_vel_x, object_vel_y, object_vel_z, object_x, object_y, object_z, object_force_x, object_force_y,
					object_force_z, object_vec[0].name, t, object_vec[0].num_nodes, object_vec[0].num_tets, object_tet_connectivity, object_nodal_area, false, object_node_curvature);
			}
			else {
				tecplot.tecplot_output_lagrangian_object_gpu(object_vel_x, object_vel_y, object_vel_z, object_x, object_y, object_z, object_force_x, object_force_y, object_force_z, globals, domain, t
					, object_vec[0].name, object_vec[0].num_nodes, object_vec[0].depth_nodes, object_vec[0].radial_nodes);
			}


			soln.clone(temp_soln);
			//only output at decreasing order of magnitudes - save space on hard drive
			/*if (convergence_residual.max_error() < pow(10, output_residual_threshold)) {*/
			tecplot.tecplot_output_unstructured_soln(globals, Mesh, soln, bcs, t, pp, residual_worker, delta_t_local, visc_factor);


			output_residual_threshold = output_residual_threshold - 1;
			soln.output(globals.output_file, globals, domain);


			square << < n_node_Blocks, blockSize >> > (total_object_nodes, object_force_x);
			square << < n_node_Blocks, blockSize >> > (total_object_nodes, object_force_y);
			square << < n_node_Blocks, blockSize >> > (total_object_nodes, object_force_z);

			total<128> << < n_node_Blocks, blockSize >> > (object_force_x, force_x_block, total_object_nodes);
			total<128> << < n_node_Blocks, blockSize >> > (object_force_y, force_y_block, total_object_nodes);
			total<128> << < n_node_Blocks, blockSize >> > (object_force_z, force_z_block, total_object_nodes);
			post_kernel_checks();
			cudaDeviceSynchronize();
			convergence_residual.force_totals(force_x_block, force_y_block, force_z_block, n_node_Blocks, total_object_nodes);
	

			error_output << t << ", " << convergence_residual.max_error() << ", " <<
				convergence_residual.rho_rms << ", " << convergence_residual.u_rms << ", " <<
				convergence_residual.v_rms << ", " <<
				convergence_residual.w_rms << " , FMG cycle: " << fmg << endl;

			cout << "time t=" << time << " error e =" << convergence_residual.max_error()
				<< " delta_t:" << delta_t <<  " Stopped: " << stopped << std::endl; 
			cout << "drag: " << drag_t1 << endl;

			cout << "time t=" << time << " force f =" << convergence_residual.force_x 
				<< " , " << convergence_residual.force_y << " , " << convergence_residual.force_z << std::endl;
			cout << "time t=" << time << " Area Coefficient =" << area_coefficient
				<< " , Volume Coefficient = " << volume_coefficient << " , MAD Coefficient" << MAD_coefficient << std::endl;
			//max_u << t << "," << soln.get_u(center_node) << "," << force << endl;
			

			cout << "max node : " << object_vec[0].max_index << endl;
			cout << "min node : " << object_vec[0].min_index << endl;

			cout << "max node velocity: " << object_vel_x[object_vec[0].max_index] << endl;
			cout << "min node velocity: " << object_vel_x[object_vec[0].min_index] << endl;
				

		}
		
		

	/*	if (convergence_residual.max_error() < local_tolerance || time > td || stopped) {*/
		if (convergence_residual.max_error() < local_tolerance || time > td || (stopped && globals.testcase == 6)) {
			if (mg == 0) {
				soln.clone(temp_soln);
				cout << "convergence" << endl;
				cout << "time t=" << time << " error e =" << convergence_residual.max_error()
					<< " delta_t:" << delta_t << " Stopped: " << stopped <<  std::endl;
				error_output.close();
				energy_output.close();
				debug_log.close();
				vortex_output.close();
				max_u.close();

				// vortex calcs
				soln.update_unstructured_bcs(bcs, Mesh, domain, t);
				grads.Get_LS_Gradients(bcs, Mesh, domain, soln, globals);
				pp.cylinder_post_processing(Mesh, globals, grads, bcs, soln, domain,res_face);
				// pp.calc_vorticity(x_gradients,y_gradients);
				  //pp.calc_streamfunction(Mesh,globals,bcs);
				tecplot.tecplot_output_unstructured_soln(globals, Mesh, soln, bcs, t, pp, residual_worker, delta_t_local, visc_factor);


				if (object_vec[0].type == 2) {
					tecplot.spring_network_output_gpu(globals.output_file, false, object_vel_x, object_vel_y, object_vel_z, object_x, object_y, object_z, object_force_x, object_force_y,
						object_force_z, object_vec[0].name, timesteps, object_vec[0].num_nodes, object_vec[0].num_tets, object_tet_connectivity, object_nodal_area,true, object_node_curvature);
					tecplot.spring_network_output_gpu(globals.output_file, false, object_vel_x, object_vel_y, object_vel_z, object_x, object_y, object_z, object_force_x, object_force_y,
						object_force_z, object_vec[0].name, timesteps, object_vec[0].num_nodes, object_vec[0].num_tets, object_tet_connectivity, object_nodal_area, false, object_node_curvature);
				}
				else {
					tecplot.tecplot_output_lagrangian_object_gpu(object_vel_x, object_vel_y, object_vel_z, object_x, object_y, object_z, object_force_x, object_force_y, object_force_z, globals, domain, timesteps
						, object_vec[0].name, object_vec[0].num_nodes, object_vec[0].depth_nodes, object_vec[0].radial_nodes);
				}

				//soln.output_centrelines(globals.output_file,globals,Mesh,time);
			}
			cudaProfilerStop();

			return;
		}


	}


	//    pp.calc_vorticity(x_gradients,y_gradients);
		//pp.calc_streamfunction(Mesh,globals,bcs);
	cudaProfilerStop();
	soln.clone(temp_soln);
	cout << "out of time" << endl;
	error_output.close();
	energy_output.close();
	vortex_output.close();
	debug_log.close();
	max_u.close();
	//pp.cylinder_post_processing(Mesh, globals, grads, bcs, soln, domain, wall_shear_stress);
	tecplot.tecplot_output_unstructured_soln(globals, Mesh, soln, bcs, timesteps, pp, residual_worker, delta_t_local, visc_factor);

	if (object_vec[0].type == 2) {
		tecplot.spring_network_output_gpu(globals.output_file, false, object_vel_x, object_vel_y, object_vel_z, object_x, object_y, object_z, object_force_x, object_force_y,
			object_force_z, object_vec[0].name, timesteps, object_vec[0].num_nodes, object_vec[0].num_tets, object_tet_connectivity, object_nodal_area,true, object_node_curvature);
		tecplot.spring_network_output_gpu(globals.output_file, false, object_vel_x, object_vel_y, object_vel_z, object_x, object_y, object_z, object_force_x, object_force_y,
			object_force_z, object_vec[0].name, timesteps, object_vec[0].num_nodes, object_vec[0].num_tets, object_tet_connectivity, object_nodal_area, false, object_node_curvature);
	}
	else {
		tecplot.tecplot_output_lagrangian_object_gpu(object_vel_x, object_vel_y, object_vel_z, object_x, object_y, object_z, object_force_x, object_force_y, object_force_z, globals, domain, timesteps
			, object_vec[0].name, object_vec[0].num_nodes, object_vec[0].depth_nodes, object_vec[0].radial_nodes);
	}



}

void gpu_solver::check_periodic_object(int blockSize, double * object_z, double * node_add_block, int n_node_Blocks, int total_object_nodes, domain_geometry &domain, 
	int n_Cell_Blocks){

	//get centroid
	double sum_add;
	sum_add = 0.0;
	total<128> << < n_node_Blocks, blockSize >> > (object_z, node_add_block, total_object_nodes);
	post_kernel_checks();
	cudaDeviceSynchronize();

	for (int nd = 0; nd < n_node_Blocks; nd++) {
		sum_add = sum_add + node_add_block[nd];
	}

	double object_centre_z;

	object_centre_z = sum_add / total_object_nodes;

	//centre 
	if (object_centre_z > (domain.origin_z + 0.75* domain.Z)) {

		add_double << < n_Cell_Blocks, blockSize >> > (total_object_nodes, object_z, (-0.5*domain.Z));

	}
	post_kernel_checks();
	cudaDeviceSynchronize();

}

void gpu_solver::get_weighted_average(gradients &grads, int i, int neighbour, double m1, double m2,
	vector_var &u, vector_var &v, vector_var &w, vector_var &rho, unstructured_mesh &mesh)
{
	double a, b, x, y, z;

	//check for boundary condition

	//use boundary cell gradients as these are at cell face
	if (neighbour > mesh.get_n_cells()) {
		x = grads.get_u(neighbour).x;
		y = grads.get_u(neighbour).y;
		z = grads.get_u(neighbour).z;
		u.set_equal(x, y, z);

		x = grads.get_v(neighbour).x;
		y = grads.get_v(neighbour).y;
		z = grads.get_v(neighbour).z;
		v.set_equal(x, y, z);


		x = grads.get_w(neighbour).x;
		y = grads.get_w(neighbour).y;
		z = grads.get_w(neighbour).z;
		w.set_equal(x, y, z);


		x = grads.get_rho(neighbour).x;
		y = grads.get_rho(neighbour).y;
		z = grads.get_rho(neighbour).z;
		rho.set_equal(x, y, z);


	}
	else {

		a = m1 + m2;
		b = m2 / a;
		a = m1 / a;


		x = grads.get_u(i).x * a + grads.get_u(neighbour).x *b;
		y = grads.get_u(i).y * a + grads.get_u(neighbour).y *b;
		z = grads.get_u(i).z * a + grads.get_u(neighbour).z *b;
		u.set_equal(x, y, z);

		x = grads.get_v(i).x * a + grads.get_v(neighbour).x *b;
		y = grads.get_v(i).y * a + grads.get_v(neighbour).y *b;
		z = grads.get_v(i).z * a + grads.get_v(neighbour).z *b;
		v.set_equal(x, y, z);


		x = grads.get_w(i).x * a + grads.get_w(neighbour).x *b;
		y = grads.get_w(i).y * a + grads.get_w(neighbour).y *b;
		z = grads.get_w(i).z * a + grads.get_w(neighbour).z *b;
		w.set_equal(x, y, z);


		x = grads.get_rho(i).x * a + grads.get_rho(neighbour).x *b;
		y = grads.get_rho(i).y * a + grads.get_rho(neighbour).y *b;
		z = grads.get_rho(i).z * a + grads.get_rho(neighbour).z *b;
		rho.set_equal(x, y, z);
	}


}


vector_var gpu_solver::get_e_alpha(int k, double &lattice_weight, double c, double PI) {

	vector_var temp;
	int x, y, z;
	//get e_alpha again
	if (k > 0 && k < 5) { //

		x = round(cos((k - 1)*PI / 2) * c);
		y = round(sin((k - 1)*PI / 2)* c);
		z = 0; //update in 3D
		lattice_weight = 1.0 / 9.0;
	}
	else if (k > 4) {

		x = round(sqrt(2) * cos((k - 5)*PI / 2 + PI / 4) * c);
		y = round(sqrt(2) * sin((k - 5)*PI / 2 + PI / 4) * c);
		z = 0; //update in 3D
		lattice_weight = 1.0 / 36.0;

	}
	else {
		x = 0;
		y = 0;
		z = 0;
		lattice_weight = 4.0 / 9.0;
	}
	temp.x = x;
	temp.y = y;
	temp.z = z;


	return temp;
}

void gpu_solver::populate_e_alpha(vector<vector_var> &e_alpha, double *lattice_weight, double c, double PI, int j) {

	vector_var temp;
	int x[15] = { 0,1,-1,0,0,0,0,1,-1, 1,-1,1,-1,-1,1 };
	int y[15] = { 0,0,0,1,-1,0,0,1,-1,1,-1,-1,1,1,-1 };
	int z[15] = { 0,0,0,0,0,1,-1,1,-1,-1,1,1,-1,1,-1 };
	//get e_alpha again

	for (int k = 0; k < j; k++) {
		if (k > 0 && k < 7) { //

			lattice_weight[k] = 1.0 / 9.0;

		}
		else if (k > 6) {


			lattice_weight[k] = 1.0 / 72.0;

		}
		else {

			lattice_weight[k] = 2.0 / 9.0;
		}



		temp.x = x[k];
		temp.y = y[k];
		temp.z = z[k];

		e_alpha.push_back(temp);


	}



}

void gpu_solver::get_cell_gradients(Mesh &Mesh, int i, int neighbour, int j, Solution &temp_soln,
	vector_var &delta_rho, vector_var &delta_rho1,
	vector_var &delta_u, vector_var &delta_u1,
	vector_var &delta_v, vector_var &delta_v1,
	Boundary_Conditions &bcs) {


	int neighbour_1, neighbour_2;
	vector_var cell_1, cell_2;
	// is it N-S or E-W
	if (j == 2) {


		neighbour_1 = Mesh.get_w_node(i);
		neighbour_2 = Mesh.get_e_node(i);

	}
	else {
		neighbour_1 = Mesh.get_s_node(i);
		neighbour_2 = Mesh.get_n_node(i);

	}

	// get neighbouring cells of cells
	Mesh.get_centroid(neighbour_1, cell_1);
	Mesh.get_centroid(neighbour_2, cell_2);

	delta_rho.Get_Gradient(temp_soln.get_rho(neighbour_1), temp_soln.get_rho(neighbour_2)
		, cell_1, cell_2);
	delta_u.Get_Gradient(temp_soln.get_u(neighbour_1), temp_soln.get_u(neighbour_2)
		, cell_1, cell_2);
	delta_v.Get_Gradient(temp_soln.get_v(neighbour_1), temp_soln.get_v(neighbour_2)
		, cell_1, cell_2);


	// get gradient of neighbouring cell
	if (j == 2) {

		neighbour_1 = Mesh.get_w_node(neighbour);
		neighbour_2 = Mesh.get_e_node(neighbour);

	}
	else {
		neighbour_1 = Mesh.get_s_node(neighbour);
		neighbour_2 = Mesh.get_n_node(neighbour);

	}

	// get neighbouring cells of cells
	Mesh.get_centroid(neighbour_1, cell_1);
	Mesh.get_centroid(neighbour_2, cell_2);

	delta_rho1.Get_Gradient(temp_soln.get_rho(neighbour_1), temp_soln.get_rho(neighbour_2)
		, cell_1, cell_2);
	delta_u1.Get_Gradient(temp_soln.get_u(neighbour_1), temp_soln.get_u(neighbour_2)
		, cell_1, cell_2);
	delta_v1.Get_Gradient(temp_soln.get_v(neighbour_1), temp_soln.get_v(neighbour_2)
		, cell_1, cell_2);

}

void gpu_solver::cell_interface_variables(int j, int i, vector_var &interface_node, int &neighbour, double &interface_area,
	vector_var &cell_normal, Boundary_Conditions &boundary_conditions, bc_var &bc,
	Mesh &Mesh, vector_var &cell_2) {

	switch (j) {

	case 0: // West
		interface_node.x = Mesh.get_west_x(i);
		interface_node.y = Mesh.get_west_y(i);
		interface_node.z = Mesh.get_west_z(i);
		neighbour = Mesh.get_w_node(i);
		interface_area = Mesh.get_w_area(i);
		cell_normal.x = Mesh.get_w_i(i);
		cell_normal.y = Mesh.get_w_j(i);
		cell_normal.z = Mesh.get_w_k(i);
		break;

	case 1: // South
		interface_node.x = Mesh.get_south_x(i);
		interface_node.y = Mesh.get_south_y(i);
		interface_node.z = Mesh.get_south_z(i);
		neighbour = Mesh.get_s_node(i);
		interface_area = Mesh.get_s_area(i);
		cell_normal.x = Mesh.get_s_i(i);
		cell_normal.y = Mesh.get_s_j(i);
		cell_normal.z = Mesh.get_s_k(i);

		break;
	case 2: // East
		interface_node.x = Mesh.get_east_x(i);
		interface_node.y = Mesh.get_east_y(i);
		interface_node.z = Mesh.get_east_z(i);
		interface_area = Mesh.get_e_area(i);
		neighbour = Mesh.get_e_node(i);
		cell_normal.x = Mesh.get_e_i(i);
		cell_normal.y = Mesh.get_e_j(i);
		cell_normal.z = Mesh.get_e_k(i);

		break;
	case 3: // North
		interface_node.x = Mesh.get_north_x(i);
		interface_node.y = Mesh.get_north_y(i);
		interface_node.z = Mesh.get_north_z(i);
		neighbour = Mesh.get_n_node(i);
		interface_area = Mesh.get_n_area(i);
		cell_normal.x = Mesh.get_n_i(i);
		cell_normal.y = Mesh.get_n_j(i);
		cell_normal.z = Mesh.get_n_k(i);

		break;
	case 4: // Front
		interface_node.x = Mesh.get_front_x(i);
		interface_node.y = Mesh.get_front_y(i);
		interface_node.z = Mesh.get_front_z(i);
		neighbour = Mesh.get_f_node(i);
		interface_area = Mesh.get_f_area(i);
		cell_normal.x = Mesh.get_f_i(i);
		cell_normal.y = Mesh.get_f_j(i);
		cell_normal.z = Mesh.get_f_k(i);

		break;
	case 5: // Back
		interface_node.x = Mesh.get_back_x(i);
		interface_node.y = Mesh.get_back_y(i);
		interface_node.z = Mesh.get_back_z(i);
		neighbour = Mesh.get_b_node(i);
		interface_area = Mesh.get_b_area(i);
		cell_normal.x = Mesh.get_b_i(i);
		cell_normal.y = Mesh.get_b_j(i);
		cell_normal.z = Mesh.get_b_k(i);
		break;


	}
	//        cell_2.x = Mesh.get_centroid_x(neighbour);
	//        cell_2.y = Mesh.get_centroid_y((neighbour));
	//        cell_2.z = Mesh.get_centroid_z(neighbour);

}



void gpu_solver::cell_interface_variables(int face, int i, vector_var &interface_node, int &neighbour, double &interface_area,
	vector_var &cell_normal, Boundary_Conditions &boundary_conditions, bc_var &bc,
	unstructured_mesh &Mesh, vector_var &cell_2, vector_var &cell_1) {



	interface_node.x = Mesh.get_face_x(face);
	interface_node.y = Mesh.get_face_y(face);
	interface_node.z = Mesh.get_face_z(face);

	neighbour = Mesh.get_mesh_neighbour(face);
	interface_area = Mesh.get_face_area(face);
	cell_normal.x = Mesh.get_face_i(face);
	cell_normal.y = Mesh.get_face_j(face);
	cell_normal.z = Mesh.get_face_k(face);


	cell_2.x = Mesh.get_centroid_x(neighbour);
	cell_2.y = Mesh.get_centroid_y((neighbour));
	cell_2.z = Mesh.get_centroid_z(neighbour);


}



void gpu_solver::get_cell_nodes(std::vector<int> &cell_nodes, Boundary_Conditions &bcs, int neighbour,
	Mesh &Mesh, int i, int j) {

	//current cell
	cell_nodes.clear();
	if (bcs.get_bc(i) || bcs.get_bc(neighbour)) {
		cell_nodes.push_back(i);
		cell_nodes.push_back(neighbour);

	}
	else if (j == 2) {
		cell_nodes.push_back(i);
		cell_nodes.push_back(Mesh.get_n_node(i));
		//cell_nodes.push_back(Mesh.get_e_node(i));
		//cell_nodes.push_back(Mesh.get_w_node(i));
		cell_nodes.push_back(Mesh.get_s_node(i));
		cell_nodes.push_back(neighbour);
		cell_nodes.push_back(Mesh.get_n_node(neighbour));
		//cell_nodes.push_back(Mesh.get_e_node(neighbour));
		//cell_nodes.push_back(Mesh.get_w_node(neighbour));
		cell_nodes.push_back(Mesh.get_s_node(neighbour));
	}
	else {
		cell_nodes.push_back(i);
		//cell_nodes.push_back(Mesh.get_n_node(i));
		cell_nodes.push_back(Mesh.get_e_node(i));
		cell_nodes.push_back(Mesh.get_w_node(i));
		// cell_nodes.push_back(Mesh.get_s_node(i));
		cell_nodes.push_back(neighbour);
		//cell_nodes.push_back(Mesh.get_n_node(neighbour));
		cell_nodes.push_back(Mesh.get_e_node(neighbour));
		cell_nodes.push_back(Mesh.get_w_node(neighbour));
		//cell_nodes.push_back(Mesh.get_s_node(neighbour));

	}
}



//get CFL numbers for inviscid and viscous matrices
// see what time stepping results
void gpu_solver::populate_cfl_areas(double3 *cfl_areas, unstructured_mesh &Mesh) {

	double area_x, area_y, area_z;
	int face;
	double3 temp;
	for (int i = 0; i < Mesh.get_n_cells(); i++) {
		area_x = 0;
		area_y = 0;
		area_z = 0;

		// time step condition as per OpenFoam calcs
		for (int f = 0; f < Mesh.gradient_faces[i].size(); f++) {
			face = Mesh.gradient_faces[i][f];

			// eigen values as per Zhaoli guo(2004) - preconditioning

			//method as per Jiri Blasek: CFD Principles and Application Determination of Max time Step

			// need to calulate correct direction of face vector

			area_x = area_x + fabs(Mesh.get_face_i(face)*Mesh.get_face_area(face));
			area_y = area_y + fabs(Mesh.get_face_j(face)*Mesh.get_face_area(face));
			area_z = area_z + fabs(Mesh.get_face_k(face)*Mesh.get_face_area(face));

		}

		temp.x = area_x / 2;
		temp.y = area_y / 2;
		temp.z = area_z / 2;
		cfl_areas[i] = temp;
	}

	return;
}




//get CFL numbers for inviscid and viscous matrices
// see what time stepping results



void gpu_solver::inverse_weighted_distance_interpolation(double &u, double &v, double &rho, Boundary_Conditions &bcs,
	Mesh &Mesh, domain_geometry &domain, Solution &soln, vector_var &interface_node,
	int k, int i, int neighbour, vector<vector_var> &e_alpha, int j, std::vector<int> &cell_nodes) {

	// get interface node
	double w_u, w_v, w_rho, w_sum, w;  // weighted macros

	// get 8 nodes'
	w_u = 0.0;
	w_v = 0.0;
	w_rho = 0.0;
	w_sum = 0.0;

	double r;
	r = 0.0;
	double dt;
	if (j == 2) {
		dt = Mesh.get_delta_t_e(i);
	}
	else {
		dt = Mesh.get_delta_t_n(i);
	}

	//get displacements
	vector_var node_displacement, target_node;

	// get target node
	target_node.x = interface_node.x - e_alpha[k].x * dt;
	target_node.y = interface_node.y - e_alpha[k].y * dt;
	target_node.z = interface_node.z - e_alpha[k].z * dt;

	// comment code as it's unused and incompatible with KAY architecture
	//for (auto &it : cell_nodes) {
	//	node_displacement.x = Mesh.get_centroid_x(it) - target_node.x;
	//	node_displacement.y = Mesh.get_centroid_y(it) - target_node.y;
	//	node_displacement.z = Mesh.get_centroid_z(it) - target_node.z;

	//	r = node_displacement.Magnitude();

	//	//
	//	if (r < 10e-5) {
	//		u = soln.get_u(it);
	//		w = soln.get_v(it);
	//		rho = soln.get_rho(it);
	//		return;

	//	}

	//	//get weight for this cc
	//	w = pow(1 / r, 2.0);

	//	// sum weighted cc values
	//	w_u = w_u + w * soln.get_u(it);
	//	w_v = w_v + w * soln.get_v(it);
	//	w_rho = w_rho + w * soln.get_rho(it);
	//	w_sum = w_sum + w;

	//}

	// calc u v rho for target node
	u = w_u / w_sum;
	v = w_v / w_sum;
	rho = w_rho / w_sum;

}

void gpu_solver::find_real_time(double* delta_t_local, double* local_time, bool* calc_face,
	unstructured_mesh &Mesh, bool* calc_cell) {

	// for each cell check cell calc check if time is greater than neighbouring cells;
	int nb;

	for (int i = 0; i < Mesh.get_total_cells(); i++) {
		// first update actual time
		if (calc_cell[i]) {
			local_time[i] = local_time[i] + delta_t_local[i];
		}
	}


	for (int i = 0; i < Mesh.get_total_cells(); i++) {


		calc_cell[i] = true;
		for (int j = 0; j < Mesh.gradient_cells[i].size(); j++) {
			nb = Mesh.gradient_cells[i][j];
			if (local_time[i] > local_time[nb]) {
				calc_cell[i] = false;
				j = Mesh.gradient_cells[i].size();
			}
		}
	}

	// then for each face calc if it should be calculated
	for (int k = 0; k < Mesh.get_n_faces(); k++) {
		calc_face[k] = false;
		if (calc_cell[Mesh.get_mesh_owner(k)] || calc_cell[Mesh.get_mesh_neighbour(k)]) {
			calc_face[k] = true;
		}
	}




}

void gpu_solver::post_kernel_checks() {

	cudaDeviceSynchronize();
	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess)
	{
		// print the CUDA error message and exit
		printf("CUDA error: %s\n", cudaGetErrorString(error));
		exit(-1);
	}


}


template <typename T>
void gpu_solver::bcs_to_array(T* target, Boundary_Conditions &bcs, int total_nodes, std::string name) {

	for (int i = 0; i < total_nodes; i++) {
		if (name.compare("vel_type") == 0) {
			target[i] = bcs.get_vel_type(i);

		}
		else if (name.compare("rho_type") == 0) {
			target[i] = bcs.get_rho_type(i);

		}
	}
}



template <typename T>
void gpu_solver::bcs_to_array_double(T* target, Boundary_Conditions &bcs, int total_nodes, std::string name) {

	for (int i = 0; i < total_nodes; i++) {
		if (name.compare("bcs") == 0) {

			double4 temp;
			temp.x = bcs.get_u(i);
			temp.y = bcs.get_v(i);
			temp.z = bcs.get_w(i);
			temp.w = bcs.get_rho(i);
			target[i] = temp;

		}

	}
}


void gpu_solver::lagrangian_object_to_array(std::vector<lagrangian_object> &obj_vec, double* &x_ref, double* &y_ref, double* &z_ref, double* &x, double* &y, double* & z,
	double* &x0, double* &y0, double* & z0, int * & tet_connectivity) {

	int n = 0;
	for (int i = 0; i < obj_vec.size(); i++) {
		for (int j = 0; j < obj_vec[i].num_nodes; j++) {

			x_ref[n] = obj_vec[i].node_x_ref[j];
			x[n] = obj_vec[i].node_x[j];
			x0[n] = obj_vec[i].node_x[j];

			y_ref[n] = obj_vec[i].node_y_ref[j];
			y[n] = obj_vec[i].node_y[j];
			y0[n] = obj_vec[i].node_y[j];

			z_ref[n] = obj_vec[i].node_z_ref[j];
			z[n] = obj_vec[i].node_z[j];
			z0[n] = obj_vec[i].node_z[j];
		

			n++;
		}

		// transfer tet related variables
		for (int k = 0; k < obj_vec[i].num_tets * 3; k++) {

			tet_connectivity[k] = obj_vec[i].tet_connectivity[k];

		}

	}
}


void gpu_solver::lagrangian_spring_network_to_array(std::vector<lagrangian_object> &obj_vec, double* &x_ref, double* &y_ref, double* &z_ref, double* &x, double* &y, double* & z,
	double* &x0, double* &y0, double* & z0, int * & tet_connectivity, int * & node_neighbours,
	int * & node_freedom, int * &  spring_connectivity , int * & tri_neighbours, double* loading_factors , int * & tri_neighbours_minor, int * & node_neighbours_minor) {

	int n = 0;
	for (int i = 0; i < obj_vec.size(); i++) {
		for (int j = 0; j < obj_vec[i].num_nodes; j++) {

			x_ref[n] = obj_vec[i].sorted_coordinate[j*3 ];
			x[n] = obj_vec[i].sorted_coordinate[j * 3];
			x0[n] = obj_vec[i].sorted_coordinate[j * 3];

			y_ref[n] = obj_vec[i].sorted_coordinate[j * 3 +1];
			y[n] = obj_vec[i].sorted_coordinate[j * 3 + 1];
			y0[n] = obj_vec[i].sorted_coordinate[j * 3 + 1];

			z_ref[n] = obj_vec[i].sorted_coordinate[j * 3 + 2];
			z[n] = obj_vec[i].sorted_coordinate[j * 3 + 2];
			z0[n] = obj_vec[i].sorted_coordinate[j * 3 + 2];
			
			loading_factors[n] = obj_vec[i].external_loading_factor[j];

			node_freedom[n] = obj_vec[i].sorted_nodFreedom[j];

			for (int k = 0; k < obj_vec[i].maxFreedom; k++) {
				node_neighbours[n* obj_vec[i].maxFreedom + k] = obj_vec[i].nodIdx[j*obj_vec[i].maxFreedom + k];
				tri_neighbours[n* obj_vec[i].maxFreedom + k] = obj_vec[i].tri_neighbours[j*obj_vec[i].maxFreedom + k];
				tri_neighbours_minor[n* obj_vec[i].maxFreedom + k] = obj_vec[i].tri_neighbours_minor[j*obj_vec[i].maxFreedom + k];
				node_neighbours_minor[n* obj_vec[i].maxFreedom + k] = obj_vec[i].minor_node_neighbour[j*obj_vec[i].maxFreedom + k];

				/*tri_neighbours[n* obj_vec[i].maxFreedom + k] = 0;
				tri_neighbours_minor[n* obj_vec[i].maxFreedom + k] = 0;
				node_neighbours_minor[n* obj_vec[i].maxFreedom + k] = 0; */
				
			}
			
			n++;
		}

		// transfer tet related variables
		for (int k = 0; k < obj_vec[i].num_tets * 3; k++) {

			tet_connectivity[k] = obj_vec[i].sorted_triIdx[k];

		}

		// transfer spring related variables
		for (int k = 0; k < obj_vec[i].num_springs * 4; k++) {

			spring_connectivity[k] = obj_vec[i].sprIdx[k];

		}

	}
}



void gpu_solver::lagrangian_spring_network_to_array_restart(std::vector<lagrangian_object> &obj_vec, double* &x_force, double* &y_force, double* &z_force, double* &x, double* &y, double* & z,
	double* &vel_x, double* &vel_y, double* & vel_z, double *&loading_factors) {

	int n = 0;
	for (int i = 0; i < obj_vec.size(); i++) {
		for (int j = 0; j < obj_vec[i].num_nodes; j++) {

			x[n] = obj_vec[i].node_x[j];
			x_force[n] = obj_vec[i].node_force_x[j];
			vel_x[n] = obj_vec[i].node_vel_x[j];

			y[n] = obj_vec[i].node_y[j];
			y_force[n] = obj_vec[i].node_force_y[j];
			vel_y[n] = obj_vec[i].node_vel_x[j];

			z[n] = obj_vec[i].node_z[j];
			z_force[n] = obj_vec[i].node_force_z[j];
			vel_z[n] = obj_vec[i].node_vel_z[j];
			loading_factors[n] = obj_vec[i].external_loading_factor[j];

			n++;
		}


	}
}


template <typename T>
void gpu_solver::mesh_to_array(T* target, unstructured_mesh &mesh, int total_nodes, std::string name, int testcase) {


	for (int i = 0; i < total_nodes; i++) {
		if (name.compare("volume") == 0) {
			target[i] = mesh.get_cell_volume(i);

		}
		else if (name.compare("gradient_stencil") == 0) {
			for (int j = 0; j < 6; j++) {
				target[i * 6 + j] = mesh.gradient_cells[i][j];
			}

		}
		else if (name.compare("mesh_owner") == 0) {
			target[i] = mesh.get_mesh_owner(i);

		}
		else if (name.compare("mesh_neighbour") == 0) {
			target[i] = mesh.get_mesh_neighbour(i);

		}
		else if (name.compare("surface_area") == 0) {
			target[i] = mesh.get_face_area(i);

		}
		else if (name.compare("streaming_dt") == 0) {
			target[i] = mesh.get_delta_t_face(i);

		}
		else if (name.compare("plateau_x") == 0) {
			if (testcase != 2) {
				target[i] = mesh.get_plateau_x(i);
			}
		}
		else if (name.compare("plateau_y") == 0) {
			if (testcase != 2) {
				target[i] = mesh.get_plateau_y(i);
			}
		}

		else if (name.compare("plateau_z") == 0) {
			if (testcase != 2) {
				target[i] = mesh.get_plateau_z(i);
			}
		}

		

	}

}



template <typename T>
void gpu_solver::mesh_to_array_double(T* target, unstructured_mesh &mesh, int total_nodes, std::string name)
{
	for (int i = 0; i < total_nodes; i++) {
		if (name.compare("cell_centroid") == 0) {

			double3 temp;
			temp.x = mesh.get_centroid_x(i);
			temp.y = mesh.get_centroid_y(i);
			temp.z = mesh.get_centroid_z(i);

			target[i] = temp;

		}
		else if (name.compare("face_normal") == 0) {

			double3 temp;
			temp.x = mesh.get_face_i(i);
			temp.y = mesh.get_face_j(i);
			temp.z = mesh.get_face_k(i);

			target[i] = temp;
		}
		else if (name.compare("face_centroid") == 0) {

			double3 temp;
			temp.x = mesh.get_face_x(i);
			temp.y = mesh.get_face_y(i);
			temp.z = mesh.get_face_z(i);

			target[i] = temp;
		}


	}
}

template <typename T>
void gpu_solver::gradients_to_array_double(T* target, gradients &grads, int total_nodes, std::string name)
{
	for (int i = 0; i < total_nodes; i++) {
		if (name.compare("RHS_array") == 0) {

			for (int j = 0; j < 6; j++) {


				double3 temp;
				temp.x = double(grads.RHS_x[i * 6 + j]);
				temp.y = double(grads.RHS_y[i * 6 + j]);
				temp.z = double(grads.RHS_z[i * 6 + j]);

				target[i * 6 + j] = temp;

			}


		}
	}
}




template <typename T>
void gpu_solver::gradients_to_array(T* target, gradients &grads, int total_nodes, std::string name)
{



	for (int i = 0; i < total_nodes; i++) {
		if (name.compare("LHS_xx") == 0) {
			target[i] = double(grads.LHS_xx[i]);
		}
		else if (name.compare("LHS_xy") == 0) {
			target[i] = grads.LHS_xy[i];
		}
		else if (name.compare("LHS_xz") == 0) {
			target[i] = grads.LHS_xz[i];
		}
		else if (name.compare("LHS_yx") == 0) {
			target[i] = grads.LHS_yx[i];
		}
		else if (name.compare("LHS_yy") == 0) {
			target[i] = grads.LHS_yy[i];
		}
		else if (name.compare("LHS_yz") == 0) {
			target[i] = grads.LHS_yz[i];
		}
		else if (name.compare("LHS_zx") == 0) {
			target[i] = grads.LHS_zx[i];
		}
		else if (name.compare("LHS_zy") == 0) {
			target[i] = grads.LHS_zy[i];
		}
		else if (name.compare("LHS_zz") == 0) {
			target[i] = grads.LHS_zz[i];
		}
	}
}






void gpu_solver::soln_to_double(double4* target, Solution &soln_a, int total_nodes) {


	for (int i = 0; i < total_nodes; i++) {
		double4 tmp;
		tmp.w = soln_a.get_rho(i);
		tmp.x = soln_a.get_u(i);
		tmp.y = soln_a.get_v(i);
		tmp.z = soln_a.get_w(i);
		target[i] = tmp;
	}


};

void gpu_solver::update_viscosities(double *object_x, double * object_y, double * object_z, double * cell_block, int n_cell_Blocks, int n_node_Blocks,  unstructured_mesh &Mesh, double * visc_factor, 
	double3 *cell_centroid, double num_nodes, double * node_block, int blockSize, double delta_h, double pi, double internal_viscosity_ratio, int num_tets, int * tet_connectivity) {

	//find min x,y,z

	double min_x, min_y, min_z;
	double max_x, max_y, max_z;

	min<128> << < n_node_Blocks, blockSize >> > (object_x, node_block,num_nodes);
	post_kernel_checks();
	cudaDeviceSynchronize();
	min_x = min_block(node_block, n_node_Blocks);

	min<128> << < n_node_Blocks, blockSize >> > (object_y, node_block, num_nodes);
	post_kernel_checks();
	cudaDeviceSynchronize();
	min_y = min_block(node_block, n_node_Blocks);

	min<128> << < n_node_Blocks, blockSize >> > (object_z, node_block, num_nodes);
	post_kernel_checks();
	cudaDeviceSynchronize();
	min_z = min_block(node_block, n_node_Blocks);

	// find max x, y, z

	max<128> << < n_node_Blocks, blockSize >> > (object_x, node_block, num_nodes);
	post_kernel_checks();
	cudaDeviceSynchronize();
	max_x = max_block(node_block, n_node_Blocks);

	max<128> << < n_node_Blocks, blockSize >> > (object_y, node_block, num_nodes);
	post_kernel_checks();
	cudaDeviceSynchronize();
	max_y = max_block(node_block, n_node_Blocks);

	max<128> << < n_node_Blocks, blockSize >> > (object_z, node_block, num_nodes);
	post_kernel_checks();
	cudaDeviceSynchronize();
	max_z = max_block(node_block, n_node_Blocks);

	fill_double << < n_cell_Blocks, blockSize >> > (Mesh.get_n_cells(), visc_factor, 1.0);

	/*update_interior_viscosities << < n_cell_Blocks, blockSize >> > (Mesh.get_n_cells(), cell_centroid, object_x, object_y, object_z,
		min_x, max_x, min_y, max_y, min_z, max_z, num_nodes, visc_factor, delta_h, pi, internal_viscosity_ratio);
*/
	post_kernel_checks();
	cudaDeviceSynchronize();

	update_interior_viscosities_plucker << < n_cell_Blocks, blockSize >> > (Mesh.get_n_cells(), cell_centroid, object_x, object_y, object_z,
		min_x, max_x, min_y, max_y, min_z, max_z, num_nodes, visc_factor, delta_h, pi, internal_viscosity_ratio, num_tets, tet_connectivity);

	post_kernel_checks();
	cudaDeviceSynchronize();

	return;
}

	


void gpu_solver::run_spring_network_functions(int total_object_nodes, double * object_force_x, double *object_force_y, double * object_force_z, double *object_x, double * object_y, double * object_z,
	std::vector<lagrangian_object> &object_vec , int * object_node_neighbours,int * object_node_freedom, double *object_tet_area, 
	double * tet_normal_x, double *  tet_normal_y, double * tet_normal_z, int * object_tri_neighbours, double * spring_contour_lengths, double *  spring_constants,
	double *object_node_curvature, double * object_node_curvature_0, double * spring_angle, double * object_vel_x, double *object_vel_y, double *object_vel_z, double * tet_add_block,
	int n_tet_Blocks, int blockSize, double * node_add_block, int n_node_Blocks, int total_object_tets, int* object_tet_connectivity,
	double *object_nodal_area, double *object_tet_volume, double *object_tet_area_0, int * object_tri_neighbours_minor, int * object_node_neighbours_minor, bool initials, double * spring_constants_pow,
	double * spring_angle_0, double * object_nodal_area_0) {


	update_spring_network_tet_parameters << < n_tet_Blocks, blockSize >> > (total_object_tets, object_x, object_y, object_z,
		object_tet_connectivity, object_nodal_area, object_tet_volume, object_tet_area, tet_normal_x, tet_normal_y, tet_normal_z);

	//sum volume in reduction
	post_kernel_checks();
	cudaDeviceSynchronize();

	total<128> << < n_tet_Blocks, blockSize >> > (object_tet_volume, tet_add_block, total_object_tets);
	post_kernel_checks();
	cudaDeviceSynchronize();

	object_vec[0].add_volume(tet_add_block, initials, n_tet_Blocks);

	//clone tet area for intitials
	if (initials) {
		clone_a_to_b << < n_tet_Blocks, blockSize >> > (total_object_tets, object_tet_area, object_tet_area_0); // area_0 holds tet area at start of simulation
		post_kernel_checks();
		cudaDeviceSynchronize();
		
	}
	
	//sum area in reduction
	total<128> << < n_tet_Blocks, blockSize >> > (object_tet_area, tet_add_block, total_object_tets);
	post_kernel_checks();
	cudaDeviceSynchronize();
	object_vec[0].add_area(tet_add_block, initials, n_tet_Blocks);
		   
	preprocess_spring_network_spring_parameters << < n_node_Blocks, blockSize >> > (total_object_nodes, object_node_neighbours, object_node_freedom, object_x, object_y, object_z, object_vec[0].maxFreedom,
		spring_contour_lengths, spring_constants, object_nodal_area, object_node_curvature, object_tri_neighbours, object_tet_area, tet_normal_x, tet_normal_y, tet_normal_z, spring_angle,
		object_vec[0].shear_modulus, object_vec[0].wlc_ratio, initials, spring_constants_pow);
	post_kernel_checks();
	cudaDeviceSynchronize();

	if (initials) {

		get_spring_contour_length << < n_node_Blocks, blockSize >> > (total_object_nodes, object_node_neighbours, object_node_freedom, object_x, object_y, object_z, object_vec[0].maxFreedom,
			spring_contour_lengths,  object_vec[0].wlc_ratio);
		post_kernel_checks();
		cudaDeviceSynchronize();
		clone_a_to_b << < n_node_Blocks, blockSize >> > (total_object_nodes, object_nodal_area, object_nodal_area_0); // area_0 holds tet area at start of simulation
		post_kernel_checks();
		cudaDeviceSynchronize();
	}

	//get MAD - node curvature hasn't yet divided by area for actual curvature

	/// MAD = sigma theta *(L)
	
	total<128> << < n_node_Blocks, blockSize >> > (object_node_curvature, node_add_block, total_object_nodes);
	post_kernel_checks();
	cudaDeviceSynchronize();

	object_vec[0].add_MAD(node_add_block, initials, n_node_Blocks);

	/// curvature = theta *(L)/ Area
	divide << < n_node_Blocks, blockSize >> > (total_object_nodes,object_node_curvature, object_nodal_area );
	post_kernel_checks();
	cudaDeviceSynchronize();

	if (initials) {
		fill_double << < n_node_Blocks, blockSize >> > (total_object_nodes, object_node_curvature_0, 0.0);
		//clone_a_to_b << < n_node_Blocks, blockSize >> > (total_object_nodes, object_node_curvature, object_node_curvature_0); //curavature_0 holds curvature at start of simulation
		post_kernel_checks();
		cudaDeviceSynchronize();

		//clone_a_to_b << < n_node_Blocks, blockSize >> > (total_object_nodes, spring_angle, spring_angle_0); //spring_angle_0 holds curvature at start of simulation
		//post_kernel_checks();
		//cudaDeviceSynchronize();
	}

	update_spring_network_forces << < n_node_Blocks, blockSize >> > (total_object_nodes, object_force_x, object_force_y, object_force_z, object_x, object_y, object_z,
		object_vec[0].volume_modulus, object_vec[0].total_volume, object_vec[0].total_volume_0, object_node_neighbours, object_node_freedom, object_vec[0].maxFreedom,
		object_vec[0].local_area_modulus, object_vec[0].global_area_modulus, object_vec[0].total_area, object_vec[0].total_area_0, object_tet_area, object_tet_area_0,
		tet_normal_x, tet_normal_y, tet_normal_z, object_tri_neighbours, spring_contour_lengths, spring_constants,
		object_vec[0].global_bending_modulus, object_vec[0].total_MAD, object_vec[0].total_MAD_0, object_vec[0].membrane_thickness, object_vec[0].local_bending_modulus, object_node_curvature,
		object_node_curvature_0, spring_angle, object_vec[0].viscosity, object_vel_x, object_vel_y, object_vel_z, object_vec[0].wlc_ratio, object_tri_neighbours_minor, object_node_neighbours_minor,
		spring_constants_pow, spring_angle_0);


	post_kernel_checks();
	cudaDeviceSynchronize();
	   
	

}




void gpu_solver::run_spring_network_functions_atomics(int total_object_nodes, double * object_force_x, double *object_force_y, double * object_force_z, double *object_x, double * object_y, double * object_z,
	std::vector<lagrangian_object> &object_vec, int * object_node_neighbours, int * object_node_freedom, double *object_tet_area,
	double * tet_normal_x, double *  tet_normal_y, double * tet_normal_z, int * object_tri_neighbours, double * spring_contour_lengths, double *  spring_constants,
	double *object_node_curvature, double * object_node_curvature_0, double * spring_angle, double * object_vel_x, double *object_vel_y, double *object_vel_z, double * tet_add_block,
	int n_tet_Blocks, int blockSize, double * node_add_block, int n_node_Blocks, int total_object_tets, int* object_tet_connectivity,
	double *object_nodal_area, double *object_tet_volume, double *object_tet_area_0, int * object_tri_neighbours_minor, int * object_node_neighbours_minor, bool initials, double * spring_constants_pow,
	double * spring_angle_0, double * object_nodal_area_0, int * spring_connectivity, int total_object_springs, int n_spring_blocks, double *spring_area, double * spring_area_0) {


	update_spring_network_tet_parameters << < n_tet_Blocks, blockSize >> > (total_object_tets, object_x, object_y, object_z,
		object_tet_connectivity, object_nodal_area, object_tet_volume, object_tet_area, tet_normal_x, tet_normal_y, tet_normal_z);

	//sum volume in reduction
	post_kernel_checks();
	cudaDeviceSynchronize();

	total<128> << < n_tet_Blocks, blockSize >> > (object_tet_volume, tet_add_block, total_object_tets);
	post_kernel_checks();
	cudaDeviceSynchronize();

	object_vec[0].add_volume(tet_add_block, initials, n_tet_Blocks);

	//sum area in reduction
	total<128> << < n_tet_Blocks, blockSize >> > (object_tet_area, tet_add_block, total_object_tets);
	post_kernel_checks();
	cudaDeviceSynchronize();
	object_vec[0].add_area(tet_add_block, initials, n_tet_Blocks);

	post_kernel_checks();
	cudaDeviceSynchronize();

	get_spring_forces_atomics << < n_spring_blocks, blockSize >> > (total_object_springs, object_force_x, object_force_y, object_force_z, object_x, object_y, object_z,
		object_vec[0].volume_modulus, object_vec[0].total_volume, object_vec[0].total_volume_0, spring_connectivity, object_node_freedom, object_vec[0].maxFreedom,
		object_vec[0].local_area_modulus, object_vec[0].global_area_modulus, object_vec[0].total_area, object_vec[0].total_area_0, spring_area, spring_area_0,
		tet_normal_x, tet_normal_y, tet_normal_z, object_tri_neighbours, spring_contour_lengths, spring_constants,
		object_vec[0].global_bending_modulus, object_vec[0].total_MAD, object_vec[0].total_MAD_0, object_vec[0].membrane_thickness, object_vec[0].local_bending_modulus, object_node_curvature,
		object_node_curvature_0, spring_angle, object_vec[0].viscosity, object_vel_x, object_vel_y, object_vel_z, object_vec[0].wlc_ratio, object_tri_neighbours_minor, object_node_neighbours_minor,
		spring_constants_pow, spring_angle_0, object_vec[0].shear_modulus,
		 initials, object_vec[0].spontaneous_angle);

	post_kernel_checks();
	cudaDeviceSynchronize();

	//get MAD
	sum_product_curvature<128> << < n_node_Blocks, blockSize >> > (object_node_curvature, spring_area ,node_add_block, total_object_nodes);
	post_kernel_checks();
	cudaDeviceSynchronize();

	object_vec[0].add_MAD(node_add_block, initials, n_node_Blocks);
	
}




//get CFL numbers for inviscid and viscous matrices
// see what time stepping results
void gpu_solver::get_cfl(double &delta_t, Solution &soln
	, unstructured_mesh &Mesh, global_variables &globals, double* delta_t_local, int* delta_t_frequency, Solution &cfl_areas) {


	double factor;


	double area_x_eigen, visc_eigen;
	factor = globals.time_marching_step;

	double visc_constant;
	visc_constant = 4;

	double min_delta_t, temp;

	double effective_speed_of_sound;
	//effective_speed_of_sound = 1/sqrt(3);
	effective_speed_of_sound = globals.max_velocity* sqrt(1 - globals.pre_conditioned_gamma + pow(globals.pre_conditioned_gamma / sqrt(3) / globals.max_velocity, 2));
	//loop through cells

	min_delta_t = 100000000000;

	for (int i = 0; i < Mesh.get_n_cells(); i++) {
		delta_t_frequency[i] = 1;

		// eigen values as per Zhaoli guo(2004) - preconditioning

		  //estimation of spectral radii s per Jiri Blasek: CFD Principles and Application Determination of Max time Step

		area_x_eigen = 0;
		area_x_eigen = (fabs(soln.get_u(i)) + effective_speed_of_sound)*cfl_areas.get_u(i)
			+ (fabs(soln.get_v(i)) + effective_speed_of_sound)*cfl_areas.get_v(i)
			+ (fabs(soln.get_w(i)) + effective_speed_of_sound)*cfl_areas.get_w(i);

		area_x_eigen = area_x_eigen / globals.pre_conditioned_gamma;


		//reducing preconditioning increases viscous flux - increases eigenvalue
		visc_eigen = 2 * globals.visc / globals.pre_conditioned_gamma / soln.get_rho(i) / Mesh.get_cell_volume(i);
		visc_eigen = visc_eigen * (cfl_areas.get_u(i)*cfl_areas.get_u(i) + cfl_areas.get_v(i)*cfl_areas.get_v(i) + cfl_areas.get_w(i)* cfl_areas.get_w(i));

		area_x_eigen = area_x_eigen + visc_constant * visc_eigen;

		// use smallest time step allowed
		temp = factor * Mesh.get_cell_volume(i) / area_x_eigen;
		if (temp < 0) {
			min_delta_t = temp;


		}
		if (temp < min_delta_t) {
			min_delta_t = temp;
		}

		if (globals.time_stepping == "local" || globals.time_stepping == "talts") {
			delta_t_local[i] = temp;

		}
		else { //constant user defined time step
			delta_t_local[i] = factor;

		}
	}

	if (globals.time_stepping == "min") {
		std::fill_n(delta_t_local, Mesh.get_n_cells(), min_delta_t);
	}


	if (globals.time_stepping == "talts") {

		for (int i = 0; i < Mesh.get_n_cells(); i++) {
			delta_t_frequency[i] = pow(2, floor(log2(delta_t_local[i] / min_delta_t)));
			delta_t_local[i] = min_delta_t * delta_t_frequency[i];
		}

	}

	return;
}


void gpu_solver::find_available_device(int deviceNumber)
{


	int num_gpus = 0; //number of gpus
	cudaGetDeviceCount(&num_gpus);
	//##ERROR handling
	if (num_gpus < 1) //check if cuda device is found
	{
		throw std::runtime_error("no CUDA capable devices detected");
	}
	else if (num_gpus < deviceNumber) //check if device can be selected by deviceNumber
	{
		std::cerr << "no CUDA device " << deviceNumber << ", only " << num_gpus << " devices found" << std::endl;
		throw std::runtime_error("CUDA capable devices can't be selected");
	}


	int maxTries = num_gpus;

	cudaDeviceProp devProp;
	cudaError rc;
	checkCudaErrors(cudaGetDeviceProperties(&devProp, deviceNumber));

	/* if the gpu compute mode is set to default we use the given `deviceNumber` */
	if (devProp.computeMode == cudaComputeModeDefault)
		maxTries = 1;

	for (int deviceOffset = 0; deviceOffset < maxTries; ++deviceOffset)
	{
		const int tryDeviceId = (deviceOffset + deviceNumber) % num_gpus;
		rc = cudaSetDevice(tryDeviceId);

		if (rc == cudaSuccess)
		{
			cudaStream_t stream;
			/* \todo: Check if this workaround is needed
			 *
			 * - since NVIDIA change something in driver cudaSetDevice never
			 * return an error if another process already use the selected
			 * device if gpu compute mode is set "process exclusive"
			 * - create a dummy stream to check if the device is already used by
			 * an other process.
			 * - cudaStreamCreate fail if gpu is already in use
			 */
			rc = cudaStreamCreate(&stream);
		}

		if (rc == cudaSuccess)
		{
			cudaDeviceProp dprop;
			checkCudaErrors(cudaGetDeviceProperties(&dprop, tryDeviceId));
			printf(
				"> Set device to %d  : %c ",tryDeviceId, dprop.name);

			
			if (cudaErrorSetOnActiveProcess == cudaSetDeviceFlags(cudaDeviceScheduleSpin))
			{
				cudaGetLastError(); //reset all errors
				/* - because of cudaStreamCreate was called cudaSetDeviceFlags crashed
				 * - to set the flags reset the device and set flags again
				 */
				checkCudaErrors(cudaDeviceReset());
				checkCudaErrors(cudaSetDeviceFlags(cudaDeviceScheduleSpin));
			}
			checkCudaErrors(cudaGetLastError());
			break;
		}
		else if (rc == cudaErrorDeviceAlreadyInUse || rc == cudaErrorDevicesUnavailable)
		{
			cudaGetLastError(); //reset all errors
			

			printf(
				"> Device %d already in use, try next.", tryDeviceId);

			continue;
		}
		else
		{
			checkCudaErrors(rc); /*error message*/
		}
	}
}
;
