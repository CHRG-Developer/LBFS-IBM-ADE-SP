#include "spring_network.hpp"
#include "common_kernels.hpp"
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "wlc.hpp"
#include <stdio.h>
// Utilities and system includes
#include <helper_cuda.h>  // helper function CUDA error checking and initialization
#include <helper_functions.h>  // helper for shared functions common to CUDA Samples



//__global__ void update_spring_network_forces(int total_nodes, double * force_x, double * force_y, double * force_z, double * x, double * y, double * z,
//	double * x_ref, double * y_ref, double *  z_ref,
//	double stiffness, double radius, double pi, int object_nodes, double * vel_x, double * vel_y, double * vel_z, double delta_t, double depth, double * object_area) {
//
//
//	int index = blockIdx.x * blockDim.x + threadIdx.x;
//	int stride = blockDim.x * gridDim.x;
//
//
//	for (int i = index; i < total_nodes; i += stride) {
//
//		if (i < total_nodes) {
//
//			// add in different methods for types of deformable objects later
//			/// assume no slip boundary in this test case
//			force_x[i] = (-stiffness * object_area[i]*(x[i] - x_ref[i]));
//			force_y[i] = (-stiffness * object_area[i] * (y[i] - y_ref[i]));
//			force_z[i] = (-stiffness * object_area[i] * (z[i] - z_ref[i]));
//			
//		}
//	}
//}





__global__ void update_spring_network_forces(int total_nodes, double * force_x, double * force_y, double * force_z, double * x, double * y, double * z,
	//volume variables
	double volume_modulus, double total_volume, double total_volume_0, int * node_neighbours, int * node_degree_of_freedom, int max_freedom,
	//area variables
	double local_area_modulus, double global_area_modulus, double total_area, double total_area_0, double * tet_area, double *tet_area_0, double * normal_x, double * normal_y, double *normal_z,
	int * node_tri_neighbours,
	//shear variables
	double  * contour_length, double * spring_constant,
	//bending variables
	double global_bending_modulus, double MAD, double MAD_0, double membrane_thickness, double local_bending_modulus, double * curvature, double * curvature_0, double * spring_angle, double viscosity,
	double * vel_x, double * vel_y, double * vel_z, double wlc_ratio, int * tri_neighbours_minor, int * node_neighbours_minor, double * pow_constant, double * spring_angle_0
) {


	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	__shared__ double volume_coefficient;

	volume_coefficient = -1* volume_modulus * (total_volume / total_volume_0 - 1.0);
	
	__shared__ double area_coefficient;

	area_coefficient = -1 * global_area_modulus * (total_area / total_area_0 - 1.0);


	__shared__ double MAD_coefficient;

	MAD_coefficient = -1*2.0*  global_bending_modulus * (MAD - MAD_0) / total_area_0;  /// factor of two in Osay code

	__shared__ double spring_force_constant;

	spring_force_constant = (0.25*wlc_ratio*wlc_ratio / ((wlc_ratio - 1.0)*(wlc_ratio - 1.0)) - 0.25 + 1.0 / wlc_ratio);
	
	for (int i = index; i < total_nodes; i += stride) {

		if (i < total_nodes) {

			

			int node_a, node_b, node_c, node_e;
			double3 a, b,c,d, e;  // a is node on current spring, b is plus 1, c is minus 1 and d is common node  // e is opposite minor spring node
			double3 derivative, total;
			total = set_equal(total, 0.0);
			int tri_neighbour, tri_neighbour_m1 , tri_neighbour_minor;
			double local_area_coefficient, force_coefficient;
			double local_bending_coefficient, local_bending_coefficient_2, local_bending_coefficient_minor;
			double stretch_ratio;
			double3 normal_1, normal_2, normal_3;  //normal_1 is for current tet n, normal_2 for prevous tet n-1 , normal_3 is for neighbour of tet n
			
			double3 area_derivative;
			double normal_norm;
			double length;

			double3 vel_a, vel_d;
			int index_0, index_1, index_m1;
			

			double3 volume_derivative = set_equal(volume_derivative, 0.0);
			double3 area_force = set_equal(area_force, 0.0);
			double3 spring_force = set_equal(spring_force, 0.0);
			double3 bending_force = set_equal(bending_force, 0.0);
			double3 area_bending_force = set_equal(area_bending_force, 0.0);
			double3 viscosity_force = set_equal(viscosity_force, 0.0);

			double3 dL_ds; //derivative of length to common node


			//loop through springs of each node 
			for (int n = 0; n < node_degree_of_freedom[i] ; n++) {

				//commmon  geoemtrical parameters
				index_0 = i * max_freedom + n ;
				index_1 = i * max_freedom + (n + 1) % node_degree_of_freedom[i];
				index_m1 = i * max_freedom + ( (n - 1) % node_degree_of_freedom[i] + node_degree_of_freedom[i] ) % node_degree_of_freedom[i];
				node_a = node_neighbours[index_0];
				node_b = node_neighbours[index_1];
				node_c = node_neighbours[index_m1];
				node_e = node_neighbours_minor[index_0];
				
				a.x = x[node_a];
				a.y = y[node_a];
				a.z = z[node_a];

				b.x = x[node_b];
				b.y = y[node_b];
				b.z = z[node_b];

				c.x = x[node_c];
				c.y = y[node_c];
				c.z = z[node_c];

				d.x = x[i]; 
				d.y = y[i];
				d.z = z[i];

				e.x = x[node_e];
				e.y = y[node_e];
				e.z = z[node_e];

				tri_neighbour = node_tri_neighbours[index_0];
				tri_neighbour_m1 = node_tri_neighbours[index_m1];
				tri_neighbour_minor = tri_neighbours_minor[index_0]; // gets neighbouring triangle of tri_neighbour

				normal_1.x = normal_x[tri_neighbour];
				normal_1.y = normal_y[tri_neighbour];
				normal_1.z = normal_z[tri_neighbour];

				normal_2.x = normal_x[tri_neighbour_m1];
				normal_2.y = normal_y[tri_neighbour_m1];
				normal_2.z = normal_z[tri_neighbour_m1];

				normal_3.x = normal_x[tri_neighbour_minor];
				normal_3.y = normal_y[tri_neighbour_minor];
				normal_3.z = normal_z[tri_neighbour_minor];

				normal_norm = sqrt(dot_product(normal_1,normal_1));

				//volume conservation force
				derivative = cross_product(a, b);
				volume_derivative = volume_derivative + derivative;
								
				// area conservation force 
				local_area_coefficient = -1 * local_area_modulus * (tet_area[tri_neighbour] / tet_area_0[tri_neighbour]  -1);
			
				area_derivative.x = normal_1.y* (z[node_b] - z[node_a]) + normal_1.z * (y[node_a] - y[node_b]);
				area_derivative.x = area_derivative.x / normal_norm;

				area_derivative.y = normal_1.x * (z[node_a] - z[node_b]) + normal_1.z	 * (x[node_b] - x[node_a]);
				area_derivative.y = area_derivative.y / normal_norm;
				
				area_derivative.z = normal_1.x * (y[node_b] - y[node_a]) + normal_1.y * (x[node_a] - x[node_b]);
				area_derivative.z = area_derivative.z / normal_norm;
				
				area_force = area_force + area_derivative * 0.5 *(area_coefficient + local_area_coefficient); // 0.5 for area_derivative

				/// bending force -> area contribution

				local_bending_coefficient = -1* local_bending_modulus * (curvature_0[i] * curvature_0[i] - curvature[i] * curvature[i]) / 12.0;  // 1/6 from equations, 0.5 for area_derivative

				local_bending_coefficient_2 = -1 * local_bending_modulus * (curvature[i] - curvature_0[i] ) *0.5;
				area_bending_force = area_bending_force + area_derivative * (local_bending_coefficient);

				// angle derivative
				derivative = get_angle_derivative(normal_1, normal_2, spring_angle[index_0], a, b, c);

				d =  d- a; // get length of spring
				length = sqrt(dot_product(d, d));
				dL_ds = d / length;

				derivative = derivative * length + dL_ds * spring_angle[index_0];

				bending_force = bending_force + derivative * (local_bending_coefficient_2 + MAD_coefficient);
				
				//viscosity
				vel_a.x = vel_x[node_a];
				vel_a.y = vel_y[node_a];
				vel_a.z = vel_z[node_a];

				vel_d.x = vel_x[i];
				vel_d.y = vel_y[i];
				vel_d.z = vel_z[i];

				
				vel_d = vel_d - vel_a;

				// YE SWE SO and Fedosov DPD viscosity
			/*	viscosity_force = viscosity_force  + vel_d * -12.0 / 13.0 / sqrt(3.0) * viscosity;

				viscosity_force = viscosity_force + dL_ds *dot_product(dL_ds, vel_d)* -4.0 / 13.0 / sqrt(3.0) * viscosity;*/
				 
				// OSAY viscosity
				viscosity_force = viscosity_force  - dL_ds * dot_product(dL_ds, vel_d) * 4 / sqrt(3.0) *viscosity;
									
				//shear force
				stretch_ratio = length / contour_length[index_0];

				force_coefficient = -1 * spring_constant[index_0] * ((0.25 / ((1.0 - stretch_ratio) * (1.0 - stretch_ratio))) - 0.25 + stretch_ratio - spring_force_constant);

				// POW spring force
	/*			force_coefficient = -1* ( spring_constant[index_0] * ( (1/((1.0 - stretch_ratio) * (1.0 - stretch_ratio)) ) - 1 + 4*stretch_ratio )
											+ pow_constant[index_0] / pow(length, 2.0));
				*/

				spring_force.x += (force_coefficient ) *dL_ds.x;  /// negative in overall equation balances with derivative of length for this node
				spring_force.y += (force_coefficient ) * dL_ds.y;
				spring_force.z += (force_coefficient ) * dL_ds.z;

				//minor derivative
				d = d + a;

				local_bending_coefficient_minor = -1 * local_bending_modulus * (curvature[node_a] + curvature[node_b] - curvature_0[node_a] - curvature_0[node_b]) *0.5; /// 0.5 from equations and 0.5 for average

				derivative = get_angle_derivative_minor(normal_1, normal_3,  a, b, d, e);

				a = a - b; // get length of spring between a and b
				length = sqrt(dot_product(a, a));
				derivative = derivative * length;

				bending_force = bending_force + derivative * (local_bending_coefficient_minor + MAD_coefficient);




			}

			volume_derivative = volume_derivative / 6.0;
			
			force_x[i] = volume_coefficient * volume_derivative.x;
			force_y[i] = volume_coefficient * volume_derivative.y;
			force_z[i] = volume_coefficient * volume_derivative.z;
			
			force_x[i] += area_force.x;
			force_y[i] += area_force.y;
			force_z[i] += area_force.z;

			force_x[i] += spring_force.x;
			force_y[i] += spring_force.y;
			force_z[i] += spring_force.z;

			force_x[i] += bending_force.x ;
			force_y[i] += bending_force.y;
			force_z[i] += bending_force.z;

			force_x[i] += viscosity_force.x;
			force_y[i] += viscosity_force.y;
			force_z[i] += viscosity_force.z;

			force_x[i] += area_bending_force.x;
			force_y[i] += area_bending_force.y;
			force_z[i] += area_bending_force.z;

		}
	}
}





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
) {


	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	__shared__ double volume_coefficient;

	volume_coefficient = -1 * volume_modulus * (total_volume / total_volume_0 - 1.0);

	__shared__ double area_coefficient;

	area_coefficient = -1 * global_area_modulus * (total_area / total_area_0 - 1.0);
	
	for (int i = index; i < total_tets; i += stride) {

		if (i < total_tets) {

			int node_a, node_b, node_c;
			double3 a, b, c;  // a is node on current spring, b is plus 1, c is minus 1 and d is common node  // e is opposite minor spring node
			double3  total, force_a, force_b,force_c;
			total = set_equal(total, 0.0);
			
			double local_area_coefficient;
		
			double3 normal_1;  //normal_1 is for current tet n, normal_2 for prevous tet n-1 , normal_3 is for neighbour of tet n

			double3 volume_derivative = set_equal(volume_derivative, 0.0);
			double3 area_force = set_equal(area_force, 0.0);
		

			//loop through springs of each node 
		
				//commmon  geoemtrical parameters
			node_a = tet_connectivity[i * 3];
			node_b = tet_connectivity[i * 3 + 1];
			node_c = tet_connectivity[i * 3 + 2];

			a.x = x[node_a];
			a.y = y[node_a];
			a.z = z[node_a];

			b.x = x[node_b];
			b.y = y[node_b];
			b.z = z[node_b];

			c.x = x[node_c];
			c.y = y[node_c];
			c.z = z[node_c];

			normal_1.x = normal_x[i];
			normal_1.y = normal_y[i];
			normal_1.z = normal_z[i];
						
			//volume conservation force
			force_a = cross_product(b, c) * volume_coefficient / 6.0;
			force_b = cross_product(c, a) * volume_coefficient / 6.0;
			force_c = cross_product(a, b) * volume_coefficient / 6.0;

			// area conservation force 
			local_area_coefficient = -1 * local_area_modulus * (tet_area[i] / tet_area_0[i] - 1);
			

			force_a = force_a + get_area_derivative(b, c, normal_1) * 0.5* (area_coefficient + local_area_coefficient);
			force_b = force_b + get_area_derivative(c, a, normal_1) * 0.5* (area_coefficient + local_area_coefficient);
			force_c = force_c + get_area_derivative(a, b, normal_1) * 0.5* (area_coefficient + local_area_coefficient);


			double val;
			val = force_a.x;
			myatomicAdd(&force_x[node_a], val);
			val = force_a.y;
			myatomicAdd(&force_y[node_a], val);
			val = force_a.z;
			myatomicAdd(&force_z[node_a], val);

			val = force_b.x;
			myatomicAdd(&force_x[node_b], val);
			val = force_b.y;
			myatomicAdd(&force_y[node_b], val);
			val = force_b.z;
			myatomicAdd(&force_z[node_b], val);

			val = force_c.x;
			myatomicAdd(&force_x[node_c], val);
			val = force_c.y;
			myatomicAdd(&force_y[node_c], val);
			val = force_c.z;
			myatomicAdd(&force_z[node_c], val);

		}
	}
}




__global__ void get_reference_curvatures(int total_springs, double * x, double * y, double * z,
	 int * spring_connectivity,  double * spring_area, double *spring_area_0, double * curvature, double * curvature_0, double * spring_angle, double * spring_angle_0, bool reference_sphere, 
	double curvature_ratio
) {


	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	
	for (int i = index; i < total_springs; i += stride) {

		if (i < total_springs) {

			int node_a, node_b, node_c, node_d;
			double3 a, b, c, d, l;  // a is node on current spring, b is plus 1, c is minus 1 and d is common node  // e is opposite minor spring node
			double3 temp;

			double3 normal_a, normal_b, normal_ba;  //normal_1 is for current tet n, normal_2 for prevous tet n-1 , normal_3 is for neighbour of tet n
			double length;
			double3 centroid_1, centroid_2;
			double angle_direction;
			double angle, cos_angle;
			double norm_1;
			double norm_2;

			double local_area_coefficient;
			double spring_curvature;
			double spring_area_total;

			//loop through springs of each node 

				//commmon  geoemtrical parameters
			node_a = spring_connectivity[i * 4];
			node_b = spring_connectivity[i * 4 + 1];
			node_c = spring_connectivity[i * 4 + 2];
			node_d = spring_connectivity[i * 4 + 3];

			a.x = x[node_a];
			a.y = y[node_a];
			a.z = z[node_a];

			b.x = x[node_b];
			b.y = y[node_b];
			b.z = z[node_b];

			c.x = x[node_c];
			c.y = y[node_c];
			c.z = z[node_c];

			d.x = x[node_d];
			d.y = y[node_d];
			d.z = z[node_d];

			l = c - a;

			length = sqrt(dot_product(l, l));

			//get normals
			normal_a = cross_product((b - a), (c - a));
			normal_b = cross_product((c - a), (d - a));

			double tet_volume[10];

			//get normas
			norm_1 = sqrt(dot_product(normal_a, normal_a));
			norm_2 = sqrt(dot_product(normal_b, normal_b));


			cos_angle = dot_product(normal_a, normal_b) / (norm_1* norm_2);

			//using fmin/fmax instead of if statements 
			cos_angle = fmin(cos_angle, 1.0);
			cos_angle = fmax(cos_angle, -1.0);

			centroid_1 = triangle_centroid(a, b, c);
			centroid_2 = triangle_centroid(a, c, d);

			centroid_2 = centroid_2 - centroid_1;

			///https://www.grasshopper3d.com/forum/topics/convex-or-concave-angle-between-faces?id=2985220%3ATopic%3A954247&page=1#comments 
			//check for concave/convex

			//curvature for convex is postive , concave curvature is negative

			// get angle_direction, if negative then concave
			normal_ba = normal_b - normal_a;
			angle_direction = dot_product(normal_ba, centroid_2);

			angle = acos(cos_angle);

			angle = ((angle_direction > 0) - (angle_direction < 0)) * angle;

			c = c - a; // get length of spring
			length = sqrt(dot_product(c, c));

			spring_area_total = (norm_1 + norm_2) / 6.0;
			
			curvature[i] = angle * length / spring_area_total;

			spring_area[i * 2] = norm_1;
			spring_area[i * 2 + 1] = norm_2;
			

			if (reference_sphere) {
				curvature_0[i] = curvature_ratio / 3.265483301;

				
				spring_angle_0[i] = angle;

				spring_area_0[i * 2] = norm_1;
				spring_area_0[i * 2 + 1] = norm_2;
			}
			else {
				curvature_0[i] = 0.0;
				spring_angle_0[i] = angle;

				spring_area_0[i * 2] = norm_1;
				spring_area_0[i * 2 + 1] = norm_2;
			}

		}
		
	}
}


__global__ void get_contour_lengths(int total_springs,  double * x, double * y, double * z, int * spring_connectivity,
	double  * contour_length, double * spring_constant, double wlc_ratio, double shear_modulus, bool pivkin, double reference_length

) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	__shared__ double spring_force_constant;

	spring_force_constant = (0.25*wlc_ratio*wlc_ratio / ((wlc_ratio - 1.0)*(wlc_ratio - 1.0)) - 0.25 + 1.0 / wlc_ratio);

	for (int i = index; i < total_springs; i += stride) {

		if (i < total_springs) {

			int node_a, node_b, node_c, node_d;
			double3 a, b, c, d, l;  // a is node on current spring, b is plus 1, c is minus 1 and d is common node  // e is opposite minor spring node
			double3 temp;
			double  force_coefficient;

			double3 normal_a, normal_b;  //normal_1 is for current tet n, normal_2 for prevous tet n-1 , normal_3 is for neighbour of tet n

			double length;

			double3 vel_a, vel_c;
			double3 dL_ds;
			double3 centroid_1, centroid_2;
			double angle_direction;
			double angle, cos_angle;
			double3 viscosity_force;
			double stretch_ratio;

			double local_area_coefficient;
			double spring_curvature;
			double spring_area_total;

			//loop through springs of each node 

				//commmon  geoemtrical parameters
			node_a = spring_connectivity[i * 4];
			node_b = spring_connectivity[i * 4 + 1];
			node_c = spring_connectivity[i * 4 + 2];
			node_d = spring_connectivity[i * 4 + 3];

			a.x = x[node_a];
			a.y = y[node_a];
			a.z = z[node_a];

			b.x = x[node_b];
			b.y = y[node_b];
			b.z = z[node_b];

			c.x = x[node_c];
			c.y = y[node_c];
			c.z = z[node_c];

			d.x = x[node_d];
			d.y = y[node_d];
			d.z = z[node_d];

			l = c - a;

			length = sqrt(dot_product(l, l));

			
			c = c - a; // get length of spring
			length = sqrt(dot_product(c, c));
			
			//shear force
			
				
				if (pivkin) {
					spring_constant[i] = get_pivkin_spring_constant(shear_modulus, wlc_ratio, reference_length);
					contour_length[i] = reference_length * wlc_ratio;
				}
				else {
					spring_constant[i] = get_spring_constant(length, shear_modulus, wlc_ratio);
					contour_length[i] = length * wlc_ratio;

				}

				/*
				contour_length[i] = 7.5 * pow(10, -8) * sqrt(23867.0) / sqrt(2562.0) * wlc_ratio;*/
			

		}
	}
}








__global__ void get_spring_forces_atomics(int total_springs, double * force_x, double * force_y, double * force_z, double * x, double * y, double * z,
	//volume variables
	double volume_modulus, double total_volume, double total_volume_0, int * spring_connectivity, int * node_degree_of_freedom, int max_freedom,
	//area variables
	double local_area_modulus, double global_area_modulus, double total_area, double total_area_0, double * spring_area, double *spring_area_0, double * normal_x, double * normal_y, double *normal_z,
	int * node_tri_neighbours,
	//shear variables
	double  * contour_length, double * spring_constant,
	//bending variables
	double global_bending_modulus, double MAD, double MAD_0, double membrane_thickness, double local_bending_modulus, double * curvature, double * curvature_0, double * spring_angle, double viscosity,
	double * vel_x, double * vel_y, double * vel_z, double wlc_ratio, int * tri_neighbours_minor, int * node_neighbours_minor, double * pow_constant, double * spring_angle_0, double shear_modulus,
	bool initials, double spontaneous_angle
) {


	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	__shared__ double spring_force_constant;

	spring_force_constant = (0.25*wlc_ratio*wlc_ratio / ((wlc_ratio - 1.0)*(wlc_ratio - 1.0)) - 0.25 + 1.0 / wlc_ratio);

	__shared__ double volume_coefficient;

	volume_coefficient = -1 * volume_modulus * (total_volume / total_volume_0 - 1.0)/6.0;

	__shared__ double area_coefficient;

	area_coefficient = -1 * global_area_modulus * (total_area / total_area_0 - 1.0);
	
	__shared__ double MAD_coefficient;

	MAD_coefficient = -1 * 2.0*  global_bending_modulus * (MAD - MAD_0) / total_area_0;  /// factor of two in Osay code
	//MAD_coefficient = 0.0;

	for (int i = index; i < total_springs; i += stride) {

		if (i < total_springs) {

			int node_a, node_b, node_c,node_d;
			double3 a, b, c, d, l;  // a is node on current spring, b is plus 1, c is minus 1 and d is common node  // e is opposite minor spring node
			double3  total, force_a, force_b, force_c,force_d;
			total = set_equal(total, 0.0);
			force_a = set_equal(force_a, 0.0);
			force_b = set_equal(force_b, 0.0);
			force_c = set_equal(force_c, 0.0);
			force_d = set_equal(force_d, 0.0);
			double3 temp;
			double  force_coefficient;

			double3 normal_a,normal_b, normal_ba;  //normal_1 is for current tet n, normal_2 for prevous tet n-1 , normal_3 is for neighbour of tet n

			double length;

			double3 vel_a, vel_c;
			double3 dL_ds;
			double3 centroid_1, centroid_2;
			double angle_direction;
			double angle,cos_angle;
			double3 viscosity_force;
			double stretch_ratio;

			double d_normal_1[27];
			double d_normal_2[27];
			double d_volume_1[9];
			double d_volume_2[9];
			double d_norm_1[9];
			double d_norm_2[9];
			double d_area[12];
			double norm_1;
			double norm_2;
			double d_angle_1[9];
			double d_angle_2[9];
			double d_angle[12];
			double d_curve[12];
			double d_volume[12];

			double local_area_coefficient;
			double spring_curvature;
			double spring_area_total;

			//loop through springs of each node 

				//commmon  geoemtrical parameters
			node_a = spring_connectivity[i * 4];
			node_b = spring_connectivity[i * 4 + 1];
			node_c = spring_connectivity[i * 4 + 2];
			node_d = spring_connectivity[i * 4 + 3];

			a.x = x[node_a];
			a.y = y[node_a];
			a.z = z[node_a];

			b.x = x[node_b];
			b.y = y[node_b];
			b.z = z[node_b];

			c.x = x[node_c];
			c.y = y[node_c];
			c.z = z[node_c];

			d.x = x[node_d];
			d.y = y[node_d];
			d.z = z[node_d];

			l = c - a;

			length = sqrt(dot_product(l, l));
			
			//get normals
			normal_a =  cross_product((b - a), (c - a));
			normal_b = cross_product((c - a), (d - a));
			double tet_volume[10];

			//get normas
			norm_1 = sqrt(dot_product(normal_a, normal_a));
			norm_2 = sqrt(dot_product(normal_b, normal_b));

			//get d_normals
			get_d_normals(a, b, c, normal_a, tet_volume, d_normal_1, 0, d_volume_1);
			get_d_normals(c, d, a, normal_b, tet_volume, d_normal_2, 0, d_volume_2);

			get_d_norms(d_norm_1, d_normal_1, norm_1, normal_a);
			get_d_norms(d_norm_2, d_normal_2, norm_2, normal_b);

			//get area
			for (int j = 0; j < 9; j++) d_area[j] = d_norm_1[j];
			for (int j = 0; j < 3; j++) d_area[j + 6] += d_norm_2[j];
			for (int j = 3; j < 6; j++) d_area[j + 6] = d_norm_2[j];
			for (int j = 6; j < 9; j++) d_area[j - 6] += d_norm_2[j];
			for (int j = 0; j < 12; j++) d_area[j] *= 0.166666666666666666666666666666666666666666; // 0.5 for changjng from norm to area // divide by 3 to account for area

			//get volume
			for (int j = 0; j < 9; j++) d_volume[j] = d_volume_1[j];
			for (int j = 0; j < 3; j++) d_volume[j + 6] += d_volume_2[j];
			for (int j = 3; j < 6; j++) d_volume[j + 6] = d_volume_2[j];
			for (int j = 6; j < 9; j++) d_volume[j - 6] += d_volume_2[j];
			for (int j = 0; j < 12; j++) d_volume[j] *= (1.0/18.0); // divide by 3 to account for area and 6 for formula

			cos_angle = dot_product(normal_a, normal_b) / (norm_1* norm_2) ;

			//using fmin/fmax instead of if statements 
			cos_angle = fmin(cos_angle, 1.0);
			cos_angle = fmax(cos_angle, -1.0);

			centroid_1 = triangle_centroid(a, b, c);
			centroid_2 = triangle_centroid(a, c, d);

			centroid_2 = centroid_2 - centroid_1;

			///https://www.grasshopper3d.com/forum/topics/convex-or-concave-angle-between-faces?id=2985220%3ATopic%3A954247&page=1#comments 
			//check for concave/convex

			//curvature for convex is postive , concave curvature is negative
			normal_ba = normal_b - normal_a;
			// get angle_direction, if positive  then concave
			angle_direction = dot_product(normal_a, centroid_2);

			angle = acos(cos_angle);

			angle = ((angle_direction < 0) -  (angle_direction > 0)) * angle;
			
			if (initials) {
				spring_angle_0[i] = angle;

				spring_area_0[i * 2] = norm_1;
				spring_area_0[i * 2 + 1] = norm_2;
			}


			double3 tri_u_hvector;
			tri_u_hvector = cross_product(normal_a / norm_1, l / length);

			for (int j = 0; j < 9; j++) {
				d_angle_1[j] = -1 / norm_1 * ( tri_u_hvector.x *d_normal_1[j * 3] + tri_u_hvector.y * d_normal_1[j * 3 + 1] + tri_u_hvector.z * d_normal_1[j * 3 + 2] );
			}

			tri_u_hvector = cross_product(l / length,normal_b / norm_2);

			for (int j = 0; j < 9; j++) {
				d_angle_2[j] = -1 / norm_2 * (tri_u_hvector.x *d_normal_2[j * 3] + tri_u_hvector.y * d_normal_2[j * 3 + 1] + tri_u_hvector.z * d_normal_2[j * 3 + 2]);
			}

			for (int j = 0; j < 9; j++) d_angle[j] = d_angle_1[j];
			for (int j = 0; j < 3; j++) d_angle[j + 6] += d_angle_2[j];
			for (int j = 3; j < 6; j++) d_angle[j + 6] = d_angle_2[j];
			for (int j = 6; j < 9; j++) d_angle[j - 6] += d_angle_2[j];
			
			c = c - a; // get length of spring
			length = sqrt(dot_product(c, c));
			dL_ds = c / length;

			for (int j = 0; j < 12; j++) {
				d_curve[j] = d_angle[j] * length;
			}

			d_curve[0] = d_curve[0] - dL_ds.x * angle;
			d_curve[1] = d_curve[1] - dL_ds.y * angle;
			d_curve[2] = d_curve[2] - dL_ds.z * angle;

			d_curve[6] = d_curve[6] +dL_ds.x * angle;
			d_curve[7] = d_curve[7] +dL_ds.y * angle;
			d_curve[8] = d_curve[8] +dL_ds.z * angle;

			///////////////////////////////////////////////////////////////////////////////////////////////
			///area force  - each spring accounts for two triangles - so do only one node per triangle
			local_area_coefficient = -1 * local_area_modulus* (norm_1 / spring_area_0[i * 2] - 1) + area_coefficient;
			
			force_b.x += local_area_coefficient * d_norm_1[3];
			force_b.y += local_area_coefficient * d_norm_1[4];
			force_b.z += local_area_coefficient * d_norm_1[5];

			local_area_coefficient = -1 * local_area_modulus* (norm_2 / spring_area_0[i * 2 +1] - 1) + area_coefficient;
			
			force_d.x += local_area_coefficient * d_norm_2[3];
			force_d.y += local_area_coefficient * d_norm_2[4];
			force_d.z += local_area_coefficient * d_norm_2[5];

			///////////////////////////////////////////////////////////////////////////////////////////////
			//volume force - only do first triangle due to spring duplication
			force_a.x += volume_coefficient * d_volume[0];
			force_a.y += volume_coefficient * d_volume[1];
			force_a.z += volume_coefficient * d_volume[2];

			force_b.x += volume_coefficient * d_volume[3];
			force_b.y += volume_coefficient * d_volume[4];
			force_b.z += volume_coefficient * d_volume[5];

			force_c.x += volume_coefficient * d_volume[6];
			force_c.y += volume_coefficient * d_volume[7];
			force_c.z += volume_coefficient * d_volume[8];

			force_d.x += volume_coefficient * d_volume[9];
			force_d.y += volume_coefficient * d_volume[10];
			force_d.z += volume_coefficient * d_volume[11];


			//viscosity
			vel_a.x = vel_x[node_a];
			vel_a.y = vel_y[node_a];
			vel_a.z = vel_z[node_a];

			vel_c.x = vel_x[node_c];
			vel_c.y = vel_y[node_c];
			vel_c.z = vel_z[node_c];


			vel_c = vel_c - vel_a;

			// OSAY viscosity
			viscosity_force = dL_ds * -1 * dot_product(dL_ds, vel_c) * 4 / sqrt(3.0) *viscosity;

			force_c = force_c + viscosity_force;
			force_a = force_a - viscosity_force;

			//curvature bending
			double local_bending_coefficient;
			spring_area_total = (norm_1 + norm_2) / 6.0;
			spring_area[i*2] = norm_1;
			spring_area[i * 2 + 1] = norm_2;
			curvature[i] = angle * length / spring_area_total;
			//if (initials) {
			//	//curvature_0[i] = 0.0;
			//	curvature_0[i] = curvature[i];
			//}
			//curvature_0 can be zero - from Simon Mendex
			local_bending_coefficient = -0.5* local_bending_modulus* (curvature_0[i] - curvature[i]) * (curvature_0[i] - curvature[i]);
			//local_bending_coefficient = 0.0;

			force_a.x += local_bending_coefficient * d_area[0];
			force_a.y += local_bending_coefficient * d_area[1];
			force_a.z += local_bending_coefficient * d_area[2];

			force_b.x += local_bending_coefficient * d_area[3];
			force_b.y += local_bending_coefficient * d_area[4];
			force_b.z += local_bending_coefficient * d_area[5];

			force_c.x += local_bending_coefficient * d_area[6];
			force_c.y += local_bending_coefficient * d_area[7];
			force_c.z += local_bending_coefficient * d_area[8];

			force_d.x += local_bending_coefficient * d_area[9];
			force_d.y += local_bending_coefficient * d_area[10];
			force_d.z += local_bending_coefficient * d_area[11];
			
			//d curvature stuff
			double curvature_coefficient;
			curvature_coefficient = -1*global_bending_modulus * (curvature[i] - curvature_0[i]) + MAD_coefficient;

			force_a.x += curvature_coefficient * d_curve[0];
			force_a.y += curvature_coefficient * d_curve[1];
			force_a.z += curvature_coefficient * d_curve[2];

			force_b.x += curvature_coefficient * d_curve[3];
			force_b.y += curvature_coefficient * d_curve[4];
			force_b.z += curvature_coefficient * d_curve[5];

			force_c.x += curvature_coefficient * d_curve[6];
			force_c.y += curvature_coefficient * d_curve[7];
			force_c.z += curvature_coefficient * d_curve[8];

			force_d.x += curvature_coefficient * d_curve[9];
			force_d.y += curvature_coefficient * d_curve[10];
			force_d.z += curvature_coefficient * d_curve[11];


			//shear force
			//if (initials) {
			//	spring_constant[i] = get_spring_constant(length, shear_modulus, wlc_ratio);
			//	contour_length[i] = length * wlc_ratio;


			//	/*get_pivkin_spring_constant(shear_modulus, wlc_ratio, 2562.0);
			//	contour_length[i] = 7.5 * pow(10, -8) * sqrt(23867.0) / sqrt(2562.0) * wlc_ratio;*/
			//	
			//	
			//		//get_spring_pow_constants(i, max_freedom,n, length, shear_modulus, wlc_ratio, spring_constant, spring_constants_pow);
			//}

			stretch_ratio = length / contour_length[i];
			force_coefficient = -1 * spring_constant[i] * ((0.25 / ((1.0 - stretch_ratio) * (1.0 - stretch_ratio))) - 0.25 + stretch_ratio - spring_force_constant);
			
			force_c = force_c + dL_ds * force_coefficient;
			force_a = force_a - dL_ds * force_coefficient;

			double val;
			
			val = force_a.x;
			myatomicAdd(&force_x[node_a], val);
			val = force_a.y;
			myatomicAdd(&force_y[node_a], val);
			val = force_a.z;
			myatomicAdd(&force_z[node_a], val);

			val = force_b.x;
			myatomicAdd(&force_x[node_b], val);
			val = force_b.y;
			myatomicAdd(&force_y[node_b], val);
			val = force_b.z;
			myatomicAdd(&force_z[node_b], val);

			val = force_c.x;
			myatomicAdd(&force_x[node_c], val);
			val = force_c.y;
			myatomicAdd(&force_y[node_c], val);
			val = force_c.z;
			myatomicAdd(&force_z[node_c], val);

			val = force_d.x;
			myatomicAdd(&force_x[node_d], val);
			val = force_d.y;
			myatomicAdd(&force_y[node_d], val);
			val = force_d.z;
			myatomicAdd(&force_z[node_d], val);


		}
	}
}


/// a is opposite spring node, b is +1 node and c is -1 node
//normal 1 is ab triangle and normal 2 is ac triangle

__device__ double3 get_angle_derivative(double3 normal_1 , double3 normal_2, double angle, double3 a, double3 b, double3 c) {

	double norm1, norm2;

	double inverse_sin_angle;
	norm1 = sqrt(dot_product(normal_1, normal_1));
	norm2 = sqrt(dot_product(normal_2, normal_2));

	
	double3 derivative;


	normal_1 = normal_1 / norm1;
	normal_2 = normal_2 / norm2;

	double3 coef_1, coef_2;

	coef_1 = normal_1 -  normal_2 * cos(angle) ;
	coef_2 = normal_2 - normal_1 * cos(angle) ;

	/// dNorm1 x (N2 - cos N1)
	double x, y, z;

	x= (coef_2.y * (b.z - a.z) + coef_2.z * (a.y - b.y)) / norm1 + (coef_1.y * (b.z - c.z) + coef_1.z * (c.y - b.y)) / norm2;
	y = (coef_2.x * (a.z - b.z) + coef_2.z * (b.x - a.x)) / norm1 + (coef_1.x * (c.z - b.z) + coef_1.z * (b.x - c.x)) / norm2;
	z = (coef_2.x * (b.y - a.y) + coef_2.y * (a.x - b.x)) / norm1 + (coef_1.x * (b.y - c.y) + coef_1.y * (c.x - b.x)) / norm2;

	derivative.x = x;
	derivative.y = y;
	derivative.z = z;

	inverse_sin_angle = 1 / (-1 * sin(angle));

	if (isinf(inverse_sin_angle) || isnan(inverse_sin_angle)) {
		inverse_sin_angle = 0.0;
	}
	
	derivative = derivative * inverse_sin_angle;

	return derivative;

}


//all indices are clockwise
//a is first node, b is spring opposite, c is final node
// normal_1 is ab, normal_2 is bc
__device__ double3 get_angle_derivative_atomic(double3 normal_1, double3 normal_2, double angle, double3 a, double3 b, double3 c) {

	double norm1, norm2;

	double inverse_sin_angle;
	norm1 = sqrt(dot_product(normal_1, normal_1));
	norm2 = sqrt(dot_product(normal_2, normal_2));


	double3 derivative;


	normal_1 = normal_1 / norm1;
	normal_2 = normal_2 / norm2;

	double3 coef_1, coef_2;

	coef_1 = normal_1 - normal_2 * cos(angle);
	coef_2 = normal_2 - normal_1 * cos(angle);

	/// dNorm1 x (N2 - cos N1)
	double x, y, z;

	x = (coef_2.y * (a.z - b.z) + coef_2.z * (b.y - a.y)) / norm1 + (coef_1.y * (b.z - c.z) + coef_1.z * (c.y - b.y)) / norm2;
	y = (coef_2.x * (b.z - a.z) + coef_2.z * (a.x - b.x)) / norm1 + (coef_1.x * (c.z - b.z) + coef_1.z * (b.x - c.x)) / norm2;
	z = (coef_2.x * (a.y - b.y) + coef_2.y * (b.x - a.x)) / norm1 + (coef_1.x * (b.y - c.y) + coef_1.y * (c.x - b.x)) / norm2;

	derivative.x = x;
	derivative.y = y;
	derivative.z = z;

	inverse_sin_angle = 1 / (-1 * sin(angle));

	if (isinf(inverse_sin_angle) || isnan(inverse_sin_angle)) {
		inverse_sin_angle = 0.0;
	}

	derivative = derivative * inverse_sin_angle;

	return derivative;

}




/// a is spring node, b is spring, node, c is force nodes, 
//normal 1 is abc triangle and normal 2 is acd triangle
__device__ double3 get_angle_derivative_minor_atomic(double3 normal_1, double3 normal_2, double angle, double3 a, double3 b, double3 c) {

	double norm1, norm2;

	double inverse_sin_angle;
	norm1 = sqrt(dot_product(normal_1, normal_1));
	norm2 = sqrt(dot_product(normal_2, normal_2));


	double3 derivative;

	normal_1 = normal_1 / norm1;
	normal_2 = normal_2 / norm2;

	double3  coef_2;

	coef_2 = normal_2 - normal_1 * cos(angle);
	double x, y, z;

	x= (coef_2.y * (b.z - a.z) + coef_2.z * (a.y - b.y)) / norm1;
	y = (coef_2.x * (a.z - b.z) + coef_2.z * (b.x - a.x)) / norm1;
	z = (coef_2.x * (b.y - a.y) + coef_2.y * (a.x - b.x)) / norm1;


	derivative.x = x;

	derivative.y = y;

	derivative.z = z;

	inverse_sin_angle = 1 / (-1 * sin(angle));

	if (isinf(inverse_sin_angle) || isnan(inverse_sin_angle)) {
		inverse_sin_angle = 0.0;
	}

	derivative = derivative * inverse_sin_angle;

	return derivative;

}







/// a and b are spring primary nodes
/// c is the node where the force is applied
/// d is opposite minor node
__device__ double3 get_angle_derivative_minor(double3 normal_1, double3 normal_3,  double3 a, double3 b, double3 c, double3 d) {

	double norm1, norm3;
	double angle;
	double inverse_sin_angle;
	norm1 = sqrt(dot_product(normal_1, normal_1));
	norm3 = sqrt(dot_product(normal_3, normal_3));

	double cos_angle;
	double angle_direction;
	double3 derivative;


	normal_1 = normal_1 / norm1;
	normal_3 = normal_3 / norm3;

	double3  coef_2;

	cos_angle = dot_product(normal_1, normal_3) / sqrt(dot_product(normal_1, normal_1) * dot_product(normal_3, normal_3));

	//using fmin/fmax instead of if statements 
	cos_angle = fmin(cos_angle, 1.0);
	cos_angle = fmax(cos_angle, -1.0);

	c = triangle_centroid(a, b, c);
	d = triangle_centroid(a, b, d);

	d = d - c;

	///https://www.grasshopper3d.com/forum/topics/convex-or-concave-angle-between-faces?id=2985220%3ATopic%3A954247&page=1#comments 
	//check for concave/convex

	//curvature for convex is postive , concave curvature is negative

	// get angle_direction, if negative then concave
	angle_direction = dot_product(normal_3, d);

	angle = acos(cos_angle);

	angle = ((angle_direction > 0) - (angle_direction < 0)) * angle;

	
	coef_2 = normal_3 - normal_1 * cos(angle);

	derivative.x = (coef_2.y * (b.z - a.z) + coef_2.z * (a.y - b.y)) / norm1;

	derivative.y = (coef_2.x * (a.z - b.z) + coef_2.z * (b.x - a.x)) / norm1;

	derivative.z = (coef_2.x * (b.y - a.y) + coef_2.y * (a.x - b.x)) / norm1;

	inverse_sin_angle = 1 / (-1 * sin(angle));

	if (isinf(inverse_sin_angle) || isnan(inverse_sin_angle)) {
		inverse_sin_angle = 0.0;
	}

	derivative = derivative * inverse_sin_angle;

	return derivative;

}


__global__ void update_spring_network_tet_parameters(int total_tet_nodes, double * x, double * y, double * z, int * tet_connectivity, double * object_area, double * tet_volume, double * tet_area,
	double * normal_x, double * normal_y, double *normal_z) {


	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	

	for (int i = index; i < total_tet_nodes; i += stride) {

		if (i < total_tet_nodes) {

			// get coordinates of 3 nodes in triangle
			double3 a, b, c;
			int node_a, node_b, node_c;
			node_a = tet_connectivity[i * 3];
			node_b = tet_connectivity[i * 3 + 1];
			node_c = tet_connectivity[i * 3 + 2];

			a.x = x[node_a];
			a.y = y[node_a];
			a.z = z[node_a];

			b.x = x[node_b];
			b.y = y[node_b];
			b.z = z[node_b];

			c.x = x[node_c];
			c.y = y[node_c];
			c.z = z[node_c];

			double3 normal; //normal to triangle
			
			normal = cross_product((b - a), (c - a));
			normal_x[i] = normal.x;
			normal_y[i] = normal.y;
			normal_z[i] = normal.z;


			double area = sqrt(dot_product(normal,normal)) *0.5;  /// area of triangle is 0.5 by norm of normal  // divide by 3 for each nodes area in future step
			tet_area[i] = area; 

			// need
			double3 centroid;

			centroid = triangle_centroid(a, b, c);

			// to calculate 
			tet_volume[i] = dot_product(normal, centroid) /6;

		}
	}
}


__global__ void update_spring_network_tet_parameters_derivatives(int total_tet_nodes, double * x, double * y, double * z, int * tet_connectivity, double * object_area, double * tet_volume, double * tet_area,
	double * normal_x, double * normal_y, double *normal_z, double* d_Volume, double * d_Normal) {


	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;


	for (int i = index; i < total_tet_nodes; i += stride) {

		if (i < total_tet_nodes) {

			// get coordinates of 3 nodes in triangle
			double3 a, b, c;
			int node_a, node_b, node_c;
			node_a = tet_connectivity[i * 3];
			node_b = tet_connectivity[i * 3 + 1];
			node_c = tet_connectivity[i * 3 + 2];

			a.x = x[node_a];
			a.y = y[node_a];
			a.z = z[node_a];

			b.x = x[node_b];
			b.y = y[node_b];
			b.z = z[node_b];

			c.x = x[node_c];
			c.y = y[node_c];
			c.z = z[node_c];

			double3 normal; //normal to triangle

			normal = cross_product((b - a), (c - a));
			normal_x[i] = normal.x;
			normal_y[i] = normal.y;
			normal_z[i] = normal.z;


			double area = sqrt(dot_product(normal, normal)) *0.5;  /// area of triangle is 0.5 by norm of normal  // divide by 3 for each nodes area in future step
			tet_area[i] = area;

			// need
			double3 centroid;

			centroid = triangle_centroid(a, b, c);

			// to calculate 
			tet_volume[i] = dot_product(normal, centroid) / 6;

			d_Normal[i * 27 + 0] = 0.0;
			d_Normal[i * 27 + 1] = c.z - b.z;
			d_Normal[i * 27 + 2] = b.y - c.y;

			d_Normal[i * 27 + 3] = b.z - c.z;
			d_Normal[i * 27 + 4] = 0.0;
			d_Normal[i * 27 + 5] = c.x - b.x;

			d_Normal[i * 27 + 6] = c.y - b.y;
			d_Normal[i * 27 + 7] = b.x - c.x;
			d_Normal[i * 27 + 8] = 0.0;

			d_Normal[i * 27 + 9] = 0.0;
			d_Normal[i * 27 + 10] = a.z - c.z;
			d_Normal[i * 27 + 11] = c.y - a.y;

			d_Normal[i * 27 + 12] = c.z - a.z;
			d_Normal[i * 27 + 13] = 0.0;
			d_Normal[i * 27 + 14] = a.x - c.x;

			d_Normal[i * 27 + 15] = a.y - c.y;
			d_Normal[i * 27 + 16] = c.x - a.x;
			d_Normal[i * 27 + 17] = 0.0;

			d_Normal[i * 27 + 18] = 0.0;
			d_Normal[i * 27 + 19] = b.z - a.z;
			d_Normal[i * 27 + 20] = a.y - b.y;

			d_Normal[i * 27 + 21] = a.z - b.z;
			d_Normal[i * 27 + 22] = 0.0;
			d_Normal[i * 27 + 23] = b.x - a.x;

			d_Normal[i * 27 + 24] = b.y - a.y;
			d_Normal[i * 27 + 25] = a.x - b.x;
			d_Normal[i * 27 + 26] = 0.0;

			//divide normal by 3 to get d_centroid x normal
			normal = normal / 3.0;

			d_Volume[i * 9 + 0] = normal.x +  d_Normal[i * 27 + 1] * centroid.y + d_Normal[i * 27 + 2]* centroid.z;
			d_Volume[i * 9 + 1] = normal.y + d_Normal[i * 27 + 3] * centroid.x  + d_Normal[i * 27 + 5] * centroid.z;
			d_Volume[i * 9 + 2] = normal.z + d_Normal[i * 27 + 6] * centroid.x + d_Normal[i * 27 + 7] * centroid.y ;

			d_Volume[i * 9 + 3] = normal.x  + d_Normal[i * 27 + 10] * centroid.y + d_Normal[i * 27 + 11] * centroid.z;
			d_Volume[i * 9 + 4] = normal.y + d_Normal[i * 27 + 12] * centroid.x  + d_Normal[i * 27 + 14] * centroid.z;
			d_Volume[i * 9 + 5] = normal.z + d_Normal[i * 27 + 15] * centroid.x + d_Normal[i * 27 + 16] * centroid.y ;

			d_Volume[i * 9 + 6] = normal.x + d_Normal[i * 27 + 19] * centroid.y + d_Normal[i * 27 + 20] * centroid.z;
			d_Volume[i * 9 + 7] = normal.y + d_Normal[i * 27 + 21] * centroid.x + d_Normal[i * 27 + 23] * centroid.z;
			d_Volume[i * 9 + 8] = normal.z + d_Normal[i * 27 + 24] * centroid.x + d_Normal[i * 27 + 25] * centroid.y;


		}
	}
}


__device__ void get_d_normals( double3 a, double3 b, double3 c, double3 normal , double * tet_volume, double *d_Normal, int i, double * d_Volume ){

			// get coordinates of 3 nodes in triangle
		
			double area = sqrt(dot_product(normal, normal)) *0.5;  /// area of triangle is 0.5 by norm of normal  // divide by 3 for each nodes area in future step
			
			// need
			double3 centroid;

			centroid = triangle_centroid(a, b, c);

			// to calculate 
			tet_volume[i] = dot_product(normal, centroid) / 6;

			d_Normal[i * 27 + 0] = 0.0;
			d_Normal[i * 27 + 1] = c.z - b.z;
			d_Normal[i * 27 + 2] = b.y - c.y;

			d_Normal[i * 27 + 3] = b.z - c.z;
			d_Normal[i * 27 + 4] = 0.0;
			d_Normal[i * 27 + 5] = c.x - b.x;

			d_Normal[i * 27 + 6] = c.y - b.y;
			d_Normal[i * 27 + 7] = b.x - c.x;
			d_Normal[i * 27 + 8] = 0.0;

			d_Normal[i * 27 + 9] = 0.0;
			d_Normal[i * 27 + 10] = a.z - c.z;
			d_Normal[i * 27 + 11] = c.y - a.y;

			d_Normal[i * 27 + 12] = c.z - a.z;
			d_Normal[i * 27 + 13] = 0.0;
			d_Normal[i * 27 + 14] = a.x - c.x;

			d_Normal[i * 27 + 15] = a.y - c.y;
			d_Normal[i * 27 + 16] = c.x - a.x;
			d_Normal[i * 27 + 17] = 0.0;

			d_Normal[i * 27 + 18] = 0.0;
			d_Normal[i * 27 + 19] = b.z - a.z;
			d_Normal[i * 27 + 20] = a.y - b.y;

			d_Normal[i * 27 + 21] = a.z - b.z;
			d_Normal[i * 27 + 22] = 0.0;
			d_Normal[i * 27 + 23] = b.x - a.x;

			d_Normal[i * 27 + 24] = b.y - a.y;
			d_Normal[i * 27 + 25] = a.x - b.x;
			d_Normal[i * 27 + 26] = 0.0;

			//divide normal by 3 to get d_centroid x normal
			normal = normal / 3.0;

			d_Volume[i * 9 + 0] = normal.x + d_Normal[i * 27 + 1] * centroid.y + d_Normal[i * 27 + 2] * centroid.z;
			d_Volume[i * 9 + 1] = normal.y + d_Normal[i * 27 + 3] * centroid.x + d_Normal[i * 27 + 5] * centroid.z;
			d_Volume[i * 9 + 2] = normal.z + d_Normal[i * 27 + 6] * centroid.x + d_Normal[i * 27 + 7] * centroid.y;

			d_Volume[i * 9 + 3] = normal.x + d_Normal[i * 27 + 10] * centroid.y + d_Normal[i * 27 + 11] * centroid.z;
			d_Volume[i * 9 + 4] = normal.y + d_Normal[i * 27 + 12] * centroid.x + d_Normal[i * 27 + 14] * centroid.z;
			d_Volume[i * 9 + 5] = normal.z + d_Normal[i * 27 + 15] * centroid.x + d_Normal[i * 27 + 16] * centroid.y;

			d_Volume[i * 9 + 6] = normal.x + d_Normal[i * 27 + 19] * centroid.y + d_Normal[i * 27 + 20] * centroid.z;
			d_Volume[i * 9 + 7] = normal.y + d_Normal[i * 27 + 21] * centroid.x + d_Normal[i * 27 + 23] * centroid.z;
			d_Volume[i * 9 + 8] = normal.z + d_Normal[i * 27 + 24] * centroid.x + d_Normal[i * 27 + 25] * centroid.y;


		}




__device__ void get_d_norms( double * d_norm, double *d_normal, double norm, double3 normal) {

	// get coordinates of 3 nodes in triangle
	normal = normal / norm;

	for (int i = 0; i < 9; i++) {
		d_norm[i] = d_normal[i * 3] * normal.x + d_normal[i * 3 + 1] * normal.y + d_normal[i * 3 + 2] * normal.z;
	}

}




__global__ void preprocess_spring_network_spring_parameters(int num_nodes, int * node_neighbours, int * node_degree_of_freedom, double * x, double * y, double * z, int max_freedom,
	double * contour_length, double * spring_constant, double * node_area, double * node_curvature, int *tri_neighbours, double * tet_area, double * normal_x, double * normal_y,double * normal_z,
	double * spring_angle, double shear_modulus, double wlc_ratio, bool initials, double * spring_constants_pow) {


	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;


	for (int i = index; i < num_nodes; i += stride) {
		
		if (i < num_nodes) {
			int node_b,node_c,node_d;
			double3 a, b,c,d;
			double length;
			node_area[i] = 0.0;
			node_curvature[i] = 0.0;
			int tri_neighbour, tri_neighbour_b;

			double3 normal_a, normal_b;

			double cos_angle,angle;
			double angle_direction;
			double spring_area =0.0;
			int previous;


			for (int n = 0; n < node_degree_of_freedom[i] ; n++) {
				previous = ((n - 1) % node_degree_of_freedom[i] + node_degree_of_freedom[i]) % node_degree_of_freedom[i];
				
				node_b = node_neighbours[i*max_freedom + n];
				
				node_c = node_neighbours[i*max_freedom + previous];
				node_d = node_neighbours[i*max_freedom + (n + 1) % node_degree_of_freedom[i]];

				tri_neighbour = tri_neighbours[i*max_freedom + n];
				a.x = x[i];
				a.y = y[i];
				a.z = z[i];

				b.x = x[node_b];
				b.y = y[node_b];
				b.z = z[node_b];

				b = b - a; // get length of spring
				length = sqrt(dot_product(b, b));
				
				if (initials) {
					spring_constant[i * max_freedom + n] = get_spring_constant(length, shear_modulus, wlc_ratio);

					//get_spring_pow_constants(i, max_freedom,n, length, shear_modulus, wlc_ratio, spring_constant, spring_constants_pow);
				}

				node_area[i] += tet_area[tri_neighbour] /3;

				//curvature calcs
				tri_neighbour_b = tri_neighbours[i*max_freedom + previous];

				//get angle between two 
				normal_a.x = normal_x[tri_neighbour];
				normal_a.y = normal_y[tri_neighbour];
				normal_a.z = normal_z[tri_neighbour];

				normal_b.x = normal_x[tri_neighbour_b];
				normal_b.y = normal_y[tri_neighbour_b];
				normal_b.z = normal_z[tri_neighbour_b];

				cos_angle = dot_product(normal_a, normal_b) / sqrt(dot_product(normal_a, normal_a) * dot_product(normal_b, normal_b));

				//using fmin/fmax instead of if statements 
				cos_angle = fmin(cos_angle, 1.0);
				cos_angle = fmax(cos_angle, -1.0);

				c.x = x[node_c];
				c.y = y[node_c];
				c.z = z[node_c];

				d.x = x[node_d];
				d.y = y[node_d];
				d.z = z[node_d];

				c = triangle_centroid(a, b, c);
				d = triangle_centroid(a, b, d);

				d = d - c;

				///https://www.grasshopper3d.com/forum/topics/convex-or-concave-angle-between-faces?id=2985220%3ATopic%3A954247&page=1#comments 
				//check for concave/convex

				//curvature for convex is postive , concave curvature is negative

				// get angle_direction, if negative then concave
				angle_direction = dot_product( normal_a , d);

				angle = acos(cos_angle);
				
				angle = ((angle_direction > 0) - (angle_direction < 0) ) * angle;


				node_curvature[i] += 0.5 * angle* length;
				spring_angle[i*max_freedom + n] = angle;

			}

			

		}
	}
}



__global__ void get_spring_contour_length(int num_nodes, int * node_neighbours, int * node_degree_of_freedom, double * x, double * y, double * z, int max_freedom,
	double * contour_length, double wlc_ratio) {


	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	
	for (int i = index; i < num_nodes; i += stride) {

		if (i < num_nodes) {
			int  node_b;
			double3 a, b;
			double length;
			

			for (int n = 0; n < node_degree_of_freedom[i]; n++) {
				node_b = node_neighbours[i*max_freedom + n];

				a.x = x[i];
				a.y = y[i];
				a.z = z[i];

				b.x = x[node_b];
				b.y = y[node_b];
				b.z = z[node_b];

				b = b - a; // get length of spring
				length = sqrt(dot_product(b, b));
				contour_length[i * max_freedom + n] = length * wlc_ratio;
	
			}

		}
	}
}



__global__ void external_loading(int num_nodes, double force, double * external_loading_factors, double *force_x, double *nodal_area, double total_area) {


	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;


	for (int i = index; i < num_nodes; i += stride) {

		if (i < num_nodes) {
			force_x[i] = force_x[i] + external_loading_factors[i] * force ;  /// approximation

		}
	}

}

__device__ double get_spring_constant(double length, double shear_modulus, double wlc_ratio) {

	
	//equation (by Mingzhu) relating shear modulus to spring constant.
	const double cons_WLC = (0.5*pow(1.0 -  1/ wlc_ratio , -3.0) + 1.0);

	//equation (by Dao&Fedosov) relating shear modulus to spring constant. //confirmed to be wrong.
	//const REAL cons_WLC = 0.75/((1-r_WLC)*(1-r_WLC))-0.75+4*r_WLC+0.5*r_WLC/pow(1-r_WLC,3.0);
	//const REAL cons_WLC = 0.75/((1-WLC_RATIO)*(1-WLC_RATIO))-0.75+4*WLC_RATIO+0.5*WLC_RATIO/pow(1-WLC_RATIO,3.0);
	//std::cout << "spring_constant: " << ( (4.0*length*SHEARMODULUS)/(SQRT3*cons_WLC) ) << std::endl;
	//getchar();


	return ((4.0*wlc_ratio*length*shear_modulus) / (sqrt(3.0)*cons_WLC));
	
}

__device__ double get_pivkin_spring_constant( double shear_modulus, double wlc_ratio, double length) {


	//equation (by Mingzhu) relating shear modulus to spring constant.
	const double cons_WLC = (0.5*pow(1.0 - 1 / wlc_ratio, -3.0) + 1.0);

	//equation (by Dao&Fedosov) relating shear modulus to spring constant. //confirmed to be wrong.
	//const REAL cons_WLC = 0.75/((1-r_WLC)*(1-r_WLC))-0.75+4*r_WLC+0.5*r_WLC/pow(1-r_WLC,3.0);
	//const REAL cons_WLC = 0.75/((1-WLC_RATIO)*(1-WLC_RATIO))-0.75+4*WLC_RATIO+0.5*WLC_RATIO/pow(1-WLC_RATIO,3.0);
	//std::cout << "spring_constant: " << ( (4.0*length*SHEARMODULUS)/(SQRT3*cons_WLC) ) << std::endl;
	//getchar();

	


	return ((4.0*wlc_ratio*length*shear_modulus) / (sqrt(3.0)*cons_WLC));

}


__device__ double3 get_area_derivative(double3 a, double3 b, double3 normal ) {
	
	// dA/dS = normal *d_normal / magnitude of normal
	
	double3 area_derivative;
	double normal_norm = sqrt(dot_product(normal, normal));
	
	area_derivative.x = normal.y *(b.z - a.z) + normal.z * (a.y - b.y);
	area_derivative.x = area_derivative.x / normal_norm;

	area_derivative.y = normal.x*(a.z - b.z) + normal.z*(b.x - a.x);
	area_derivative.y = area_derivative.y / normal_norm;

	area_derivative.z = normal.x*(b.y - a.y) + normal.y*(a.x - b.x);
	area_derivative.z = area_derivative.z / normal_norm;

	return area_derivative;
}
__device__ void get_spring_pow_constants(int i, int max_freedom, int n, double length, double shear_modulus, double wlc_ratio, double * spring_constant, double * pow_constant){

	//two equations
	// shear modulus equality
	//spring constant = Kb T / (4 *persistence length )

	double x;
	x= 1 / wlc_ratio;
	double alpha, beta, gamma;
	alpha = sqrt(3.0) / length;
	beta = x * 0.5 / pow(1.0 - x,3.0)  - 0.25 / pow(1.0 - x, 2.0) + 0.25;
	gamma = sqrt(3.0) * 3.0 * 0.25 / pow(length, 3.0);  // kp coefficient
	
	//force equal zero at equilibirum
	double beta2;
	beta2 = 1.0/ pow(1.0 - x, 2.0) - 1.0 + 4.0*x;

	double spring_con , pow_k, temp;

	//get the underline part first to avoid divisions
	temp = alpha * beta + gamma * pow(length, 2.0) *beta2;
	spring_con = shear_modulus / temp;

	pow_k = -1.0 * pow(length, 2.0) * spring_con * beta2;
	spring_constant[i*max_freedom + n] = spring_con;
	pow_constant[i*max_freedom + n] = pow_k;

}




__global__ void get_area_energy(int total_object_springs, double *object_spring_area, double *object_spring_area_0, double *area_energy, double global_area_modulus)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;


	for (int i = index; i < total_object_springs; i += stride) {
		if (i < total_object_springs) {
			area_energy[i] = global_area_modulus*  0.5*( (object_spring_area[i*2] - object_spring_area_0[i * 2]) * (object_spring_area[i * 2] - object_spring_area_0[i * 2]) / object_spring_area_0[i * 2]+
				(object_spring_area[i * 2+1] - object_spring_area_0[i * 2 + 1]) * (object_spring_area[i * 2 + 1] - object_spring_area_0[i * 2 + 1]) / object_spring_area_0[i * 2 + 1]);
		}
	}

}


__global__ void  get_bending_energy(int total_object_nodes, double * object_node_curvature, double * object_node_curvature_0, double* bending_energy, double * object_nodal_area_0,
	double local_bending_modulus) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;


	for (int i = index; i < total_object_nodes; i += stride) {
		if (i < total_object_nodes) {
			bending_energy[i] = 0.5* local_bending_modulus*(object_node_curvature[i] - object_node_curvature_0[i]) * (object_node_curvature[i] - object_node_curvature_0[i]) / object_nodal_area_0[i];
		}
	}
}





__global__ void get_spring_energy(int num_nodes, int * node_neighbours, int * spring_connectivity, double * x, double * y, double * z, int max_freedom,
	double * contour_length, double wlc_ratio, double * spring_energy, double * spring_constant) {


	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	__shared__ double energy_constant_1;
	energy_constant_1 = (0.25 / (1 - wlc_ratio) / (1 - wlc_ratio) - 0.25 + wlc_ratio); 

	__shared__ double energy_constant_2;
	energy_constant_2 = wlc_ratio * wlc_ratio* (3.0 - 2.0*wlc_ratio) / (1.0 - wlc_ratio) - 4 * energy_constant_1 * wlc_ratio;


	for (int i = index; i < num_nodes; i += stride) {

		if (i < num_nodes) {
			int  node_a,node_c;
			double3 a, c;
			double length;
			double stretch_ratio;
			spring_energy[i] = 0.0;
			node_a = spring_connectivity[i * 4];
			node_c = spring_connectivity[i * 4 + 2];
			
			a.x = x[node_a];
			a.y = y[node_a];
			a.z = z[node_a];

		
			c.x = x[node_c];
			c.y = y[node_c];
			c.z = z[node_c];


			c = c - a; // get length of spring
			length = sqrt(dot_product(c, c));
			stretch_ratio = length / contour_length[i ];

			spring_energy[i] = spring_energy[i] + 0.5* spring_constant[i] * fabs(0.25*contour_length[i] *
				(stretch_ratio*stretch_ratio* (3.0 - 2.0*stretch_ratio) / (1.0 - stretch_ratio) - 4 * stretch_ratio* energy_constant_1 - energy_constant_2));
				
			 //first 0.5 is to account for double counting of springs
			

		}
	}
}
