#include "immersed_boundary_method.hpp"

#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
// Utilities and system includes
#include <helper_cuda.h>  // helper function CUDA error checking and initialization
#include <helper_functions.h>  // helper for shared functions common to CUDA Samples
#include "common_kernels.hpp"



__global__ void interpolate_velocities_on_nodes(int total_nodes, double * vel_x, double * vel_y, double * vel_z, double * x, double * y, double * z,
	double3 mesh_origin, double3 mesh_lengths, double delta_x, double4* soln, int total_cells) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	__shared__ int x_cells;
	x_cells = mesh_lengths.x / delta_x;
	__shared__ int y_cells;
	y_cells = mesh_lengths.y / delta_x;
	__shared__ double inverse_delta_x_3;
	inverse_delta_x_3 = 1 / delta_x / delta_x / delta_x;
	//loop through lagrangian nodes
	for (int n = index; n < total_nodes; n += stride) {

		if (n < total_nodes) {

			vel_x[n] = 0.0;
			vel_y[n] = 0.0;
			vel_z[n] = 0.0;

			double x_ref, y_ref, z_ref;
			double4 cell_soln;
			int x_min, x_max, y_min, y_max, z_min, z_max;

			int t = 0;
			int cell_index;
			/// assume 2 point stencil for now

			//round  x_ref to nearest vertice

			x_ref = (x[n] - mesh_origin.x) / delta_x;

			x_min = (int)round(x_ref) - 1;
			x_max = x_min + 1;

			y_ref = fabs((y[n] - mesh_origin.y) / delta_x); // ANSYS mesh sorting

			y_min = (int)round(y_ref) - 1;
			y_max = y_min + 1;

			z_ref = fabs((z[n] - mesh_origin.z) / delta_x); // ANSYS mesh sorting

			z_min = (int)round(z_ref) - 1;
			z_max = z_min + 1;

			///check assumption that each cell is equal to 1
			/// otherwise need to account for this
			//weighting must add up to 1

			for (int i = x_min; i <= x_max; i++) {
				for (int j = y_min; j <= y_max; j++) {
					for (int k = z_min; k <= z_max; k++) {
						cell_index = i + (j)* x_cells + k * y_cells*x_cells;
						if (cell_index < total_cells && cell_index >= 0) {

							double dist_x = (i + 0.5 - x_ref); //distance between cell_index and node reference
							double dist_y = (j + 0.5 - y_ref); // using ANSYS Mesh sorting have to reverse Y origin
							double dist_z = (k + 0.5 - z_ref);

							double weight_x = 1 - abs(dist_x);
							double weight_y = 1 - abs(dist_y);
							double weight_z = 1 - abs(dist_z);

							cell_soln = soln[cell_index];

							/*vel_x[n] += cell_soln.x * weight_x*weight_y*weight_z*inverse_delta_x_3;
							vel_y[n] += cell_soln.y * weight_x*weight_y*weight_z*inverse_delta_x_3;
							vel_z[n] += cell_soln.z* weight_x*weight_y*weight_z*inverse_delta_x_3;*/

							vel_x[n] += cell_soln.x * weight_x*weight_y*weight_z;
							vel_y[n] += cell_soln.y * weight_x*weight_y*weight_z;
							vel_z[n] += cell_soln.z* weight_x*weight_y*weight_z;
						}
					}
				}
			}

			t = t;

		}
	}

}


__global__ void interpolate_velocities_on_nodes_cos_kernel(int total_nodes, double * vel_x, double * vel_y, double * vel_z, double * x, double * y, double * z,
	double3 mesh_origin, double3 mesh_lengths, double delta_x, double4* soln, double pi) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	__shared__ int x_cells;
	x_cells = mesh_lengths.x / delta_x;
	__shared__ double inverse_delta_x_3;
	inverse_delta_x_3 = 1 / delta_x / delta_x/delta_x;
	//loop through lagrangian nodes
	for (int n = index; n < total_nodes; n += stride) {

		if (n < total_nodes) {

			vel_x[n] = 0.0;
			vel_y[n] = 0.0;
			vel_z[n] = 0.0;

			double x_ref, y_ref;
			double4 cell_soln;
			int x_min, x_max, y_min, y_max;

			int t = 0;
			int cell_index;
			/// assume 2 point stencil for now

			//round  x_ref to nearest vertice

			x_ref = (x[n] - mesh_origin.x) / delta_x;

			x_min = (int)round(x_ref) - 2;
			x_max = x_min + 3;

			y_ref = fabs((y[n] - mesh_origin.y) / delta_x); // ANSYS mesh sorting

			y_min = (int)round(y_ref) - 2;
			y_max = y_min + 3;

			///check assumption that each cell is equal to 1
			/// otherwise need to account for this
			//weighting must add up to 1

			for (int i = x_min; i <= x_max; i++) {
				for (int j = y_min; j <= y_max; j++) {
					cell_index = i + (j)* x_cells;
					double dist_x = (i + 0.5 - x_ref); //distance between cell_index and node reference
					double dist_y = (j + 0.5 - y_ref); // using ANSYS Mesh sorting have to reverse Y origin

					double weight_x = 0.25*(1 + cos(pi* dist_x / 2));
					double weight_y = 0.25*(1 + cos(pi* dist_y / 2));
					double weight_z = 1;

					cell_soln = soln[cell_index];

					vel_x[n] += cell_soln.x * weight_x*weight_y*weight_z*inverse_delta_x_3;
					vel_y[n] += cell_soln.y * weight_x*weight_y*weight_z*inverse_delta_x_3;
					vel_z[n] += cell_soln.z* weight_x*weight_y*weight_z*inverse_delta_x_3;

				}
			}

		}
	}

}


__global__ void update_node_forces(int total_nodes, double * force_x, double * force_y, double * force_z, double * x, double * y, double * z,
	double * x_ref, double * y_ref, double *  z_ref,
	double stiffness, double radius, double pi, int object_nodes, double * vel_x, double * vel_y, double * vel_z, double delta_t, double depth) {


	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	__shared__ double area;
	area = 2 * pi * radius / object_nodes * depth;

	for (int i = index; i < total_nodes; i += stride) {

		if (i < total_nodes) {

			// add in different methods for types of deformable objects later
			/// assume no slip boundary in this test case
			force_x[i] = (-stiffness * area*(x[i] - x_ref[i]));
			force_y[i] = (-stiffness * area* (y[i] - y_ref[i]));
			force_z[i] = (-stiffness * area* (z[i] - z_ref[i]));


			////direct forcing 
			/*force_x[i] = -1*vel_x[i]  /delta_t;
			force_y[i] = -1 * vel_y[i]/delta_t ;
			force_z[i] = 0.0;*/


		}
	}
}


__global__ void update_node_positions_tweezers(int total_nodes, double * vel_x, double * vel_y, double * vel_z, double * x, double * y, double * z,
	double delta_t, int num_nodes, double point_mass, double *force_x, double *force_y, double *force_z) {


	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	__shared__ double centre_x;
	__shared__ double centre_y;


	for (int i = index; i < total_nodes; i += stride) {
		

		if (i < total_nodes) {

			vel_x[i] += force_x[i] * delta_t / point_mass;
			vel_y[i] += force_y[i] * delta_t / point_mass;
			vel_z[i] += force_z[i] * delta_t / point_mass;
			
			// add in different methods for types of deformable objects later
			//assume constant timestep
			x[i] += vel_x[i] * delta_t ;
			y[i] += vel_y[i] * delta_t ;
			z[i] += vel_z[i] * delta_t ;

		}
	}


	return;
}


__global__ void update_node_positions(int total_nodes, double * vel_x, double * vel_y, double * vel_z, double * x, double * y, double * z,
	double delta_t, int num_nodes, int rk) {


	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	__shared__ double centre_x;
	__shared__ double centre_y;


	for (int i = index; i < total_nodes; i += stride) {
		double alpha[4] = { 1.0 , 0.5, 0.5,1.0 };

		if (i < total_nodes) {
			

			// add in different methods for types of deformable objects later
			//assume constant timestep
			x[i] += vel_x[i] * delta_t ;
			y[i] += vel_y[i] * delta_t ;
			z[i] += vel_z[i] * delta_t ;

			

		}
	}


	return;
}



/// delta_h is the 
__global__ void update_interior_viscosities(int total_cells,  double3 *cell_centroid, double * object_x, double * object_y, double * object_z,
	double min_x, double max_x, double min_y, double max_y, double min_z, double max_z , int num_nodes, double * visc_factor, double delta_h, double pi, double internal_viscosity_ratio) {


	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	__shared__ double centre_x;
	__shared__ double centre_y;

	__shared__ double x_low, x_high, y_low, y_high, z_low, z_high;
	__shared__ double centroid_x, centroid_y, centroid_z;
	x_low = min_x - 2 * delta_h;
	x_high = max_x + 2 * delta_h;

	y_low = min_y - 2 * delta_h;
	y_high = max_y + 2 * delta_h;
	
	z_low = min_z - 2 * delta_h;
	z_high = max_z + 2 * delta_h;

	centroid_x = (x_low + x_high) *0.5;
	centroid_y = (y_low + y_high) *0.5;
	centroid_z = (z_low + z_high) *0.5;

	for (int i = index; i < total_cells; i += stride) {

		if (i < total_cells &&  i >0  ) {

			//find min_distance to object node
			if (cell_centroid[i].x >= x_low && cell_centroid[i].x <= x_high
				&& cell_centroid[i].y >= y_low && cell_centroid[i].y <= y_high
				&& cell_centroid[i].z >= z_low && cell_centroid[i].z <= z_high) {

				int min_node;
				double distance, min_distance;
				min_distance = 1000000000;
				for (int j = 0; j < num_nodes; j++) {
					distance = sqrt( (cell_centroid[i].x - object_x[j]) *  (cell_centroid[i].x - object_x[j]) +
						(cell_centroid[i].y - object_y[j])*(cell_centroid[i].y - object_y[j]) +
						(cell_centroid[i].z - object_z[j])*(cell_centroid[i].z - object_z[j]) );

					min_node = (distance < min_distance) * j + (distance > min_distance) * min_node;  //use booleans to avoid diversion
					min_distance = fmin(distance, min_distance);

				}

				double centroid_distance;
				centroid_distance = sqrt((cell_centroid[i].x - centroid_x)*(cell_centroid[i].x - centroid_x) +
					(cell_centroid[i].y - centroid_y)*(cell_centroid[i].y - centroid_y) +
					(cell_centroid[i].z - centroid_z)*(cell_centroid[i].z - centroid_z));

				double node_distance;
				node_distance = sqrt((object_x[min_node] - centroid_x)*(object_x[min_node] - centroid_x) +
					(object_y[min_node] - centroid_y)*(object_y[min_node] - centroid_y) +
					(object_z[min_node] - centroid_z)*(object_z[min_node] - centroid_z));

				// assign viscosity factor
				//booleans used to implement heaviside function
				min_distance = ((node_distance > centroid_distance) - (centroid_distance > node_distance)) * min_distance;



				double visc;

				
				visc = 1.0 + (internal_viscosity_ratio- 1.0) *  (0.5*(1 + min_distance * 0.5 / delta_h + 1.0 / pi * sin(pi * min_distance *0.5 / delta_h)));

				if (internal_viscosity_ratio < 1.0) {
					visc = fmin(visc, 1.0);
					visc = fmax(visc, internal_viscosity_ratio);
				}
				else {
					visc = fmax(visc, 1.0);
					visc = fmin(visc, internal_viscosity_ratio);
				}
				
				visc_factor[i] = visc;
			}
			
			
		}
	}


	return;
}

__device__ void plucker_coordinates_cell(double3 cell, double* plucker_cell, double x_boundary) {

	plucker_cell[0] = cell.x *cell.y - x_boundary * cell.y;
	plucker_cell[1] = cell.x *cell.z - x_boundary * cell.z;
	plucker_cell[2] = cell.x - x_boundary;

	//plucker 3,4,5 are equal to zero as y,z coordinates are equal

}

__device__ void plucker_coordinates_edge(double3 a, double3 b, double* plucker_edge) {

	//technically below are plucker coordinates 3,4,5 but we don't need 0,1,2 due to zero nature of cell pluckers
	plucker_edge[0] = a.y* b.z - b.y*a.z;
	plucker_edge[1] = a.z - b.z;
	plucker_edge[2] = b.y - a.y;


}

__device__ double side_operator(double* cell, double* a) {

	double side;
	//technically below are plucker coordinates 3,4,5 but we don't need 0,1,2 due to zero nature of cell pluckers
	side = cell[0] * a[1] + cell[1] * a[2] + cell[2] * a[0];

	//plucker 3,4,5 are equal to zero as y,z coordinates are equal
	return side;
}


/// delta_h is the 
__global__ void update_interior_viscosities_plucker(int total_cells, double3 *cell_centroid, double * object_x, double * object_y, double * object_z,
	double min_x, double max_x, double min_y, double max_y, double min_z, double max_z, int num_nodes, double * visc_factor, double delta_h, double pi, double internal_viscosity_ratio,
	int num_tets, int * tet_connectivity) {


	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	__shared__ double centre_x;
	__shared__ double centre_y;

	__shared__ double x_low, x_high, y_low, y_high, z_low, z_high;
	__shared__ double centroid_x, centroid_y, centroid_z;
	x_low = min_x - 2 * delta_h;
	x_high = max_x + 2 * delta_h;

	y_low = min_y - 2 * delta_h;
	y_high = max_y + 2 * delta_h;

	z_low = min_z - 2 * delta_h;
	z_high = max_z + 2 * delta_h;

	centroid_x = (x_low + x_high) *0.5;
	centroid_y = (y_low + y_high) *0.5;
	centroid_z = (z_low + z_high) *0.5;

	for (int i = index; i < total_cells; i += stride) {

		if (i < total_cells &&  i >0) {

			//find min_distance to object node
			if (cell_centroid[i].x >= x_low && cell_centroid[i].x <= x_high
				&& cell_centroid[i].y >= y_low && cell_centroid[i].y <= y_high
				&& cell_centroid[i].z >= z_low && cell_centroid[i].z <= z_high) {

				///find nearest side, in x direction from cell centroid
				
				double x_boundary = x_high;
			

				double plucker_cell[3];
				double plucker_a[3];
				double plucker_b[3];
				double plucker_c[3];
				double side_a, side_b, side_c;
				int intersection_counter = 0;
				//make line and get plucker coordinates
				plucker_coordinates_cell(cell_centroid[i], plucker_cell, x_boundary);

				int node_a, node_b, node_c;
				double3 a, b, c;

				// get distance to nearest nodes
				double distance, min_distance;
				min_distance = 1000000000;

				//loop through triangles
				for (int tet = 0; tet < num_tets; tet++) {
					//check if less/greater than cell centroid
						//commmon  geoemtrical parameters
					node_a = tet_connectivity[tet * 3];
					node_b = tet_connectivity[tet * 3 + 1];
					node_c = tet_connectivity[tet * 3 + 2];

					a.x = object_x[node_a];
					a.y = object_y[node_a];
					a.z = object_z[node_a];

					b.x = object_x[node_b];
					b.y = object_y[node_b];
					b.z = object_z[node_b];

					c.x = object_x[node_c];
					c.y = object_y[node_c];
					c.z = object_z[node_c];

					if (a.x > cell_centroid[i].x && b.x > cell_centroid[i].x && c.x > cell_centroid[i].x) {
						// get plucker coordinates of edges
						plucker_coordinates_edge(a, b, plucker_a);
						plucker_coordinates_edge(b, c, plucker_b);
						plucker_coordinates_edge(c, a, plucker_c);

						// find sideoperator
						// https://members.loria.fr/SLazard/ARC-Visi3D/Pant-project/files/plucker.html#Side
						side_a = side_operator(plucker_cell, plucker_a);
						side_b = side_operator(plucker_cell, plucker_b);
						side_c = side_operator(plucker_cell, plucker_c);


						// inside/outside check
						if ((side_a < 0 && side_b < 0 && side_c < 0) ||
							(side_a > 0 && side_b > 0 && side_c > 0)) {
							intersection_counter++;
						}

						distance = sqrt(dot_product(cell_centroid[i] - a, cell_centroid[i] - a));
						min_distance = fmin(distance, min_distance);

						distance = sqrt(dot_product(cell_centroid[i] - b, cell_centroid[i] - b));
						min_distance = fmin(distance, min_distance);

						distance = sqrt(dot_product(cell_centroid[i] - c, cell_centroid[i] - c));
						min_distance = fmin(distance, min_distance);

					}
					

				}

				bool inside = true;
				if (intersection_counter % 2 == 0) {
					inside = false;
				}

				double visc;
				double factor;

				factor = 2.0 *fmin(0.5, (min_distance * 0.5 / delta_h + 1.0 / pi * sin(pi * min_distance / delta_h) *0.5));

				if (inside) {
					visc = (internal_viscosity_ratio + 1.0) + (internal_viscosity_ratio - 1.0) *  factor;
				}
				else {
					visc = (internal_viscosity_ratio + 1.0) + (1.0 - internal_viscosity_ratio) *  factor;

				}


				

				visc_factor[i] = visc * 0.5;
			}

		}
	}


	return;
}





__global__ void update_node_positions_rk4(int total_nodes, double * vel_x, double * vel_y, double * vel_z, double * x, double * y, double * z,
	double delta_t, int num_nodes, double * x0, double * y0, double * z0) {


	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	__shared__ double centre_x;
	__shared__ double centre_y;


	for (int i = index; i < total_nodes; i += stride) {

		if (i < total_nodes &&  i >0) {
			/*if (threadIdx.x == 0) {

				centre_x = 0.0;
				centre_y = 0.0;
			}

			__syncthreads();*/

			// add in different methods for types of deformable objects later
			//assume constant timestep
			x[i] = x0[i] + vel_x[i] * delta_t;
			y[i] = y0[i] + vel_y[i] * delta_t;
			z[i] = z0[i] + vel_z[i] * delta_t;

			/*x[i] += 0.0;
			y[i] += 0.0;
			z[i] += 0.0;*/
			//reset xo
			x0[i] = x[i];
			y0[i] = y[i];
			z0[i] = z[i];

			/*	myatomicAdd(&centre_x, (x[i] / num_nodes));
				myatomicAdd(&centre_y, (y[i] / num_nodes));

			__syncthreads();
			if (threadIdx.x == 0) {


				printf("RK4 x: %4.2f %4.2f \n", centre_x, centre_y);
			}*/

		}
	}


	return;
}


// Returns index of  element closest to target in arr[] 
__device__ int find_index_closest(double arr[], int n, double target)
{
	// Corner cases 
	if (target <= arr[0])
		return 0;
	if (target >= arr[n - 1])
		return n - 1;

	// Doing binary search 
	int i = 0, j = n, mid = 0;
	while (i < j) {
		mid = (i + j) / 2;

		if (arr[mid] == target)
			return mid;

		/* If target is less than array element,
			then search in left */
		if (target < arr[mid]) {

			// If target is greater than previous 
			// to mid, return closest of two 
			if (mid > 0 && target > arr[mid - 1])
				if (target - arr[mid - 1] >= arr[mid] - target)
					return mid;
				else
					return mid - 1;

			/* Repeat for left half */
			j = mid;
		}

		// If target is greater than mid 
		else {
			if (mid < n - 1 && target < arr[mid + 1])

				if (target - arr[mid] >= arr[mid + 1] - target)
					return mid + 1;
				else
					return mid;

			// update i 
			i = mid + 1;
		}
	}

	// Only single element left after search 
	return mid;
}



__global__ void interpolate_velocities_on_nodes_non_uniform_grid(int total_nodes, double * vel_x, double * vel_y, double * vel_z, double * x, double * y, double * z,
	double3 mesh_origin, double3 mesh_lengths, double delta_x, double4* soln, int total_cells,  double PI, double3 * centroid, double * cell_volume,
	double * plateau_x, double * plateau_y, double * plateau_z) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	__shared__ int x_cells;
	x_cells = mesh_lengths.x / delta_x;
	__shared__ int y_cells;
	y_cells = mesh_lengths.y / delta_x;
	__shared__ int z_cells;
	z_cells = mesh_lengths.z / delta_x;


	
	//loop through lagrangian nodes
	for (int n = index; n < total_nodes; n += stride) {
		double  ref_index;
		if (n < total_nodes) {
			double4 cell_soln;

			vel_x[n] = 0.0;
			vel_y[n] = 0.0;
			vel_z[n] = 0.0;

			double x_ref, y_ref, z_ref;
			int x_min, x_max, y_min, y_max, z_min, z_max;

			double val;
			int t = 0;
			int cell_index;
			/// assume 2 point stencil for now

			//find relative distance to origin

					//find relative distance to origin

			x_ref = (x[n]);
			ref_index = find_index_closest(plateau_x, x_cells + 1, x_ref);

			//4 point stencil - minus 1 to convert vertex index to cell index and get mid cell of stencil - another minus 1 to get min_cell_index
			x_min = ref_index - 2;
			x_max = x_min + 3;


			y_ref = y[n]; // ANSYS mesh sorting
			ref_index = find_index_closest(plateau_y, y_cells + 1, y_ref);

			//4 point stencil - minus 1 to convert vertex index to cell index and get mid cell of stencil - another minus 1 to get min_cell_index
			y_min = ref_index - 2;
			y_max = y_min + 3;


			z_ref = z[n]; // ANSYS mesh sorting
			ref_index = find_index_closest(plateau_z, z_cells + 1, z_ref);

			//4 point stencil - minus 1 to convert vertex index to cell index and get mid cell of stencil - another minus 1 to get min_cell_index
			z_min = ref_index - 2;
			z_max = z_min + 3;

			//vandermonde matrix for non-uniform weights see Jang 2017

			//only need first column of inverse due to nature of moment calcs

			double x1, x2, x3, x4;
			x1 = centroid[x_min].x - x[n];
			x2 = centroid[x_min + 1].x - x[n];
			x3 = centroid[x_min + 2].x - x[n];
			x4 = centroid[x_max].x - x[n];

			double wx[4];
			vandemonde_inverse(x1, x2, x3, x4, wx);

			double y1, y2, y3, y4;
			y1 = centroid[y_min*x_cells].y - y[n];
			y2 = centroid[(y_min + 1)*x_cells].y - y[n];
			y3 = centroid[(y_min + 2)*x_cells].y - y[n];
			y4 = centroid[(y_max)*x_cells].y - y[n];

			double wy[4];
			vandemonde_inverse(y1, y2, y3, y4, wy);

			double z1, z2, z3, z4;
			z1 = centroid[z_min*y_cells*x_cells].z - z[n];
			z2 = centroid[(z_min + 1)* y_cells * x_cells].z - z[n];
			z3 = centroid[(z_min + 2)* y_cells*x_cells].z - z[n];
			z4 = centroid[z_max*y_cells*x_cells].z - z[n];

			double wz[4];
			vandemonde_inverse(z1, z2, z3, z4, wz);

			double weight;

			///check assumption that each cell is equal to 1
			/// otherwise need to account for this
			//weighting must add up to 1

			for (int i = x_min; i <= x_max; i++) {
				for (int j = y_min; j <= y_max; j++) {
					for (int k = z_min; k <= z_max; k++) {
						cell_index = i + (j)* x_cells + k * y_cells*x_cells;
						if (cell_index < total_cells && cell_index >= 0) {
							weight = wx[i - x_min] * wy[j - y_min] * wz[k - z_min] ;
							
							cell_soln = soln[cell_index];

							vel_x[n] += cell_soln.x * weight;
							vel_y[n] += cell_soln.y * weight;
							vel_z[n] += cell_soln.z* weight;
						}
					}
				}
			}

		}
	}

}


__global__ void spread_forces_on_non_uniform_grid(int total_nodes, double * node_force_x, double * node_force_y, double * node_force_z, double * x, double * y, double * z,
	double3 mesh_origin, double3 mesh_lengths, double delta_x,
	double * force_x, double * force_y, double * force_z, int total_cells,double PI, double3 * centroid, double * cell_volume,
	double * plateau_x, double * plateau_y, double * plateau_z){

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	__shared__ int x_cells;
	x_cells = mesh_lengths.x / delta_x;
	__shared__ int y_cells;
	y_cells = mesh_lengths.y / delta_x;
	__shared__ int z_cells;
	z_cells = mesh_lengths.z / delta_x;

	
	//loop through lagrangian nodes
	for (int n = index; n < total_nodes; n += stride) {
		double  ref_index;
		if (n < total_nodes) {


			double x_ref, y_ref, z_ref;
			int x_min, x_max, y_min, y_max, z_min, z_max;

			double val,weight;
			int t = 0;
			int cell_index;
			/// assume 2 point stencil for now

			//find relative distance to origin

			x_ref = (x[n]);
			ref_index = find_index_closest(plateau_x, x_cells + 1, x_ref);

			//4 point stencil - minus 1 to convert vertex index to cell index and get mid cell of stencil - another minus 1 to get min_cell_index
			x_min = ref_index - 2;
			x_max = x_min + 3;


			y_ref = y[n]; // ANSYS mesh sorting
			ref_index = find_index_closest(plateau_y, y_cells + 1, y_ref);

				//4 point stencil - minus 1 to convert vertex index to cell index and get mid cell of stencil - another minus 1 to get min_cell_index
			y_min = ref_index - 2;
			y_max = y_min + 3;


			z_ref = z[n]; // ANSYS mesh sorting
			ref_index = find_index_closest(plateau_z, z_cells + 1, z_ref);

			//4 point stencil - minus 1 to convert vertex index to cell index and get mid cell of stencil - another minus 1 to get min_cell_index
			z_min = ref_index - 2;
			z_max = z_min + 3;

			//vandermonde matrix for non-uniform weights see Jang 2017

			//only need first column of inverse due to nature of moment calcs

			double x1, x2, x3,x4;
			x1 = centroid[x_min].x - x[n];
			x2 = centroid[x_min + 1].x - x[n];
			x3 = centroid[x_min+ 2].x - x[n];
			x4 = centroid[x_max].x - x[n];

			double wx[4];
			vandemonde_inverse(x1, x2, x3, x4, wx);

			double y1, y2, y3, y4;
			y1 = centroid[y_min*x_cells].y - y[n];
			y2 = centroid[(y_min + 1)*x_cells].y - y[n];
			y3 = centroid[(y_min + 2)*x_cells].y - y[n];
			y4 = centroid[(y_max)*x_cells].y - y[n];

			double wy[4];
			vandemonde_inverse(y1, y2, y3, y4, wy);

			double z1, z2, z3, z4;
			z1 = centroid[z_min*y_cells*x_cells].z - z[n];
			z2 = centroid[(z_min + 1)* y_cells * x_cells].z - z[n];
			z3 = centroid[(z_min + 2)* y_cells*x_cells].z - z[n];
			z4 = centroid[z_max*y_cells*x_cells].z - z[n];

			double wz[4];
			vandemonde_inverse(z1, z2, z3, z4, wz);
			
			for (int i = x_min; i <= x_max; i++) {
				for (int j = y_min; j <= y_max; j++) {
					for (int k = z_min; k <= z_max; k++) {

						cell_index = i + (j)* x_cells + k * y_cells*x_cells;
						if (cell_index < total_cells && cell_index >= 0) {

							weight = wx[i - x_min] * wy[j - y_min] * wz[k - z_min] / cell_volume[cell_index];

							val = node_force_x[n] * weight;
							myatomicAdd(&force_x[cell_index], val);
							val = node_force_y[n] * weight;
							myatomicAdd(&force_y[cell_index], val);
							val = node_force_z[n] * weight;
							myatomicAdd(&force_z[cell_index], val);

							t++;
						}
					}
				}
			}

		}
	}

}

__device__ double reference_sin_index(double x_ref, double mesh_length, int x_cells,double PI){
		//using sqrt sin mesh
		//using sqrt sin mesh
			double fx;
			double ref_index;

			if ((x_ref / mesh_length) > 0.5) {
				fx = ((x_ref) / mesh_length - 1)* 2.0;
				fx = fx * fx;
				ref_index = x_cells - asin(fx)*x_cells / PI;

			}
			else {
			   fx = (x_ref) / mesh_length * 2.0;
			   fx = fx * fx;
			   ref_index = asin(fx)*x_cells / PI;
			}

			return ref_index;
	
	}


//only need the first column of the inverse for the weights

__device__ void vandemonde_inverse(double x, double y, double z, double weights[]) {
	
		double sub_xy, sub_xz, sub_yz;
		sub_xy = 1 / (x - y);
		sub_xz = 1 / (x - z);
		sub_yz = 1 / (y - z);
	
		
	
		// need first column of inverse to calculate weighting function
	
		weights[0] = y *z * sub_xy * sub_xz;
	
		weights[1] =  x *z *-1 * sub_xy * sub_yz;	
	
		weights[2] = x*y  * sub_xz * sub_yz;
	
	}


__device__ void vandemonde_inverse(double x, double y, double z, double w, double weights[]) {

	double sub_xy, sub_xz, sub_xw,  sub_yz, sub_yw, sub_zw ;
	sub_xy = 1 / (x - y);
	sub_xz = 1 / (x - z);
	sub_xw = 1 / (x - w);
	sub_yz = 1 / (y - z);
	sub_yw = 1 / (y - w);
	sub_zw = 1 / (z - w);



	// need first column of inverse to calculate weighting function

	weights[0] = -1* y * z *w * sub_xy * sub_xz * sub_xw;

	weights[1] =  x * z *w * sub_xy * sub_yz *sub_yw; // two minus' cancel

	weights[2] =  -1* x * y  * w * sub_xz * sub_yz * sub_zw;

	weights[3] =  x * y  * z * sub_xw * sub_yw * sub_zw;

}





__global__ void spread_forces_on_structured_grid(int total_nodes, double * node_force_x, double * node_force_y, double * node_force_z, double * x, double * y, double * z,
	double3 mesh_origin, double3 mesh_lengths, double delta_x,
	double * force_x, double * force_y, double * force_z, int total_cells) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	__shared__ int x_cells;
	x_cells = mesh_lengths.x / delta_x;
	__shared__ int y_cells;
	y_cells = mesh_lengths.y / delta_x;
	__shared__ double inverse_delta_x_3;
	inverse_delta_x_3 = 1 / delta_x / delta_x / delta_x;
	//loop through lagrangian nodes
	for (int n = index; n < total_nodes; n += stride) {

		if (n < total_nodes) {


			double x_ref, y_ref, z_ref;
			int x_min, x_max, y_min, y_max, z_min, z_max;

			double val;
			int t = 0;
			int cell_index;
			/// assume 2 point stencil for now

			//round  x_ref to nearest vertice

			x_ref = (x[n] - mesh_origin.x) / delta_x;

			x_min = (int)round(x_ref) - 1;
			x_max = x_min + 1;

			y_ref = fabs((y[n] - mesh_origin.y) / delta_x); // ANSYS mesh sorting

			y_min = (int)round(y_ref) - 1;
			y_max = y_min + 1;

			z_ref = fabs((z[n] - mesh_origin.z) / delta_x); // ANSYS mesh sorting

			z_min = (int)round(z_ref) - 1;
			z_max = z_min + 1;


			for (int i = x_min; i <= x_max; i++) {
				for (int j = y_min; j <= y_max; j++) {
					for (int k = z_min; k <= z_max; k++) {

						cell_index = i + (j)* x_cells + k * y_cells*x_cells;
						if (cell_index < total_cells && cell_index >= 0) {


							double dist_x = (i + 0.5 - x_ref); //distance between cell_index and node reference
							double dist_y = (j + 0.5 - y_ref); // using ANSYS Mesh sorting have to reverse Y origin
							double dist_z = (k + 0.5 - z_ref);

							double weight_x = 1 - abs(dist_x);
							double weight_y = 1 - abs(dist_y);
							double weight_z = 1 - abs(dist_z);
							
							val = node_force_x[n] * weight_x*weight_y*weight_z*inverse_delta_x_3;
							myatomicAdd(&force_x[cell_index], val);
							val = node_force_y[n] * weight_x*weight_y*weight_z*inverse_delta_x_3;
							myatomicAdd(&force_y[cell_index], val);
							val = node_force_z[n] * weight_x*weight_y*weight_z*inverse_delta_x_3;
							myatomicAdd(&force_z[cell_index], val);

							t++;
						}
					}
				}
			}

		}
	}

}



__global__ void spread_forces_on_structured_grid_cos_kernel(int total_nodes, double * node_force_x, double * node_force_y, double * node_force_z, double * x, double * y, double * z,
	double3 mesh_origin, double3 mesh_lengths, double delta_x,
	double * force_x, double * force_y, double * force_z, int total_cells, double pi) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	__shared__ int x_cells;
	x_cells = mesh_lengths.x / delta_x;
	__shared__ double inverse_delta_x_2;
	inverse_delta_x_2 = 1 / delta_x / delta_x;
	//loop through lagrangian nodes
	for (int n = index; n < total_nodes; n += stride) {

		if (n < total_nodes) {

			double x_ref, y_ref;
			int x_min, x_max, y_min, y_max;
			double val;
			int t = 0;
			int cell_index;
			/// assume 2 point stencil for now

			//round  x_ref to nearest vertice

			x_ref = (x[n] - mesh_origin.x) / delta_x;

			x_min = (int)round(x_ref) - 2;
			x_max = x_min + 3;

			y_ref = fabs((y[n] - mesh_origin.y) / delta_x); // ANSYS mesh sorting

			y_min = (int)round(y_ref) - 2;
			y_max = y_min + 3;

			for (int i = x_min; i <= x_max; i++) {
				for (int j = y_min; j <= y_max; j++) {
					cell_index = i + (j)* x_cells;
					if (cell_index < total_cells) {


						double dist_x = (i + 0.5 - x_ref); //distance between cell_index and node reference
						double dist_y = (j + 0.5 - y_ref); // using ANSYS Mesh sorting have to reverse Y origin

						double weight_x = 0.25*(1 + cos(pi* dist_x / 2));
						double weight_y = 0.25*(1 + cos(pi* dist_y / 2));
						double weight_z = 1;
						val = node_force_x[n] * weight_x*weight_y*weight_z*inverse_delta_x_2;
						myatomicAdd(&force_x[cell_index], val);
						val = node_force_y[n] * weight_x*weight_y*weight_z*inverse_delta_x_2;
						myatomicAdd(&force_y[cell_index], val);
						val = node_force_z[n] * weight_x*weight_y*weight_z*inverse_delta_x_2;
						myatomicAdd(&force_z[cell_index], val);

						t++;
					}

				}
			}

		}
	}

}

