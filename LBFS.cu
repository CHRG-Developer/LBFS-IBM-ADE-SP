#include "common_kernels.hpp"

#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "global_variables.h"
// Utilities and system includes
#include <helper_cuda.h>  // helper function CUDA error checking and initialization
#include <helper_functions.h>  // helper for shared functions common to CUDA Samples
#include "LBFS.hpp"



__device__ void populate_lattice_macros_uniform(double u_lattice[], double v_lattice[],
	double w_lattice[], double rho_lattice[], double3 cell_1, double3 cell,
	double3 interface_node, int i, int neighbour,
	double3 grad_u_1, double3 grad_u, double3 grad_v_1, double3  grad_v, double3  grad_w_1, double3  grad_w, double3  grad_rho_1, double3  grad_rho,
	double4 owner_soln, double4 neighbour_soln, int n_cells, double3 cell_normal, double dt, int bcs_rho_type[], int bcs_vel_type[], double4 bcs_macros[]) {

	double3 temp1, temp2; 

	///   case 0: // center node

	temp1 = cell_1 - interface_node;
	temp2 = cell - interface_node;

	double rho_i, rho_nb, u_i, u_nb, v_i, v_nb, w_i, w_nb;
	
	int nb;
	nb = neighbour - n_cells;

	if (neighbour > n_cells) {
		double4 bcs = bcs_macros[nb];
		// vboundary = v_i + grad_boundary * distance _ib

		if (bcs_rho_type[nb] == 1) {
			rho_lattice[0] = bcs.w;

		}
		else {
			rho_lattice[0] = owner_soln.w - dot_product(temp1, grad_rho);
		}

		if (bcs_vel_type[nb] == 1) {
			u_lattice[0] = bcs.x;
			v_lattice[0] = bcs.y;
			w_lattice[0] = bcs.z;

		}
		else {
			u_lattice[0] = owner_soln.x - dot_product(temp1, grad_u);
			v_lattice[0] = owner_soln.y - dot_product(temp1, grad_v);
			w_lattice[0] = owner_soln.z - dot_product(temp1, grad_w);
		}

		


	}
	else {

		rho_i = owner_soln.w - dot_product(temp1, grad_rho_1);
		rho_nb = neighbour_soln.w - dot_product(temp2, grad_rho);
		rho_lattice[0] = (rho_i + rho_nb)*0.5;

		u_i = owner_soln.x - dot_product(temp1, grad_u_1);
		u_nb = neighbour_soln.x - dot_product(temp2, grad_u);
		u_lattice[0] = (u_i + u_nb)*0.5;

		v_i = owner_soln.y - dot_product(temp1, grad_v_1);
		v_nb = neighbour_soln.y - dot_product(temp2, grad_v);
		v_lattice[0] = (v_i + v_nb)*0.5;

		w_i = owner_soln.z - dot_product(temp1, grad_w_1);
		w_nb = neighbour_soln.z - dot_product(temp2, grad_w);
		w_lattice[0] = (w_i + w_nb)*0.5;





		if ((u_lattice[0] * cell_normal.x + v_lattice[0] * cell_normal.y + w_lattice[0] * cell_normal.z) > 0.01) {


			grad_u = grad_u_1;
			grad_v = grad_v_1;
			grad_w = grad_w_1;
			grad_rho = grad_rho_1;
		}
		else if ((u_lattice[0] * cell_normal.x + v_lattice[0] * cell_normal.y + w_lattice[0] * cell_normal.z) < -0.01) {
			grad_u = grad_u;
			grad_v = grad_v;
			grad_w = grad_w;
			grad_rho = grad_rho;
		}
		else {
			grad_u.x = (grad_u_1.x + grad_u.x)*0.5;
			grad_u.y = (grad_u_1.y + grad_u.y)*0.5;
			grad_u.z = (grad_u_1.z + grad_u.z)*0.5;

			grad_v.x = (grad_v_1.x + grad_v.x)*0.5;
			grad_v.y = (grad_v_1.y + grad_v.y)*0.5;
			grad_v.z = (grad_v_1.z + grad_v.z)*0.5;

			grad_w.x = (grad_w_1.x + grad_w.x)*0.5;
			grad_w.y = (grad_w_1.y + grad_w.y)*0.5;
			grad_w.z = (grad_w_1.z + grad_w.z)*0.5;

			grad_rho.x = (grad_rho_1.x + grad_rho.x)*0.5;
			grad_rho.y = (grad_rho_1.y + grad_rho.y)*0.5;
			grad_rho.z = (grad_rho_1.z + grad_rho.z)*0.5;

		}



	}

	//set i and nb to cell interface values

	rho_i = rho_lattice[0];
	rho_nb = rho_lattice[0];

	u_i = u_lattice[0];
	u_nb = u_lattice[0];

	v_i = v_lattice[0];
	v_nb = v_lattice[0];

	w_i = w_lattice[0];
	w_nb = w_lattice[0];
	//merge gradients





	///  case 1:west_node


	rho_lattice[1] = rho_nb - grad_rho.x* dt;
	u_lattice[1] = u_nb - grad_u.x* dt;
	v_lattice[1] = v_nb - grad_v.x *dt;
	w_lattice[1] = w_nb - grad_w.x *dt;



	///  case 2: // east_node

	rho_lattice[2] = rho_nb + grad_rho.x* dt;
	u_lattice[2] = u_nb + grad_u.x* dt;
	v_lattice[2] = v_nb + grad_v.x *dt;
	w_lattice[2] = w_nb + grad_w.x *dt;



	///   case 3: // bottom node

	rho_lattice[3] = rho_nb - grad_rho.y* dt;
	u_lattice[3] = u_nb - grad_u.y* dt;
	v_lattice[3] = v_nb - grad_v.y *dt;
	w_lattice[3] = w_nb - grad_w.y *dt;



	///   case 4: // top node


	rho_lattice[4] = rho_nb + grad_rho.y* dt;
	u_lattice[4] = u_nb + grad_u.y* dt;
	v_lattice[4] = v_nb + grad_v.y *dt;
	w_lattice[4] = w_nb + grad_w.y *dt;

	///   case 5: // back node

	rho_lattice[5] = rho_nb - grad_rho.z* dt;
	u_lattice[5] = u_nb - grad_u.z* dt;
	v_lattice[5] = v_nb - grad_v.z *dt;
	w_lattice[5] = w_nb - grad_w.z *dt;

	///   case 6: // front node

	rho_lattice[6] = rho_nb + grad_rho.z* dt;
	u_lattice[6] = u_nb + grad_u.z* dt;
	v_lattice[6] = v_nb + grad_v.z *dt;
	w_lattice[6] = w_nb + grad_w.z *dt;



	/// case 7: back bottom west

	rho_lattice[7] = rho_nb - grad_rho.x* dt
		- grad_rho.y* dt
		- grad_rho.z* dt;
	u_lattice[7] = u_nb - grad_u.x* dt
		- grad_u.y* dt
		- grad_u.z* dt;
	v_lattice[7] = v_nb - grad_v.x* dt
		- grad_v.y* dt
		- grad_v.z* dt;
	w_lattice[7] = w_nb - grad_w.x* dt
		- grad_w.y* dt
		- grad_w.z* dt;

	/// case 9: front bottom west

	rho_lattice[9] = rho_nb - grad_rho.x* dt
		- grad_rho.y* dt
		+ grad_rho.z* dt;
	u_lattice[9] = u_nb - grad_u.x* dt
		- grad_u.y* dt
		+ grad_u.z* dt;
	v_lattice[9] = v_nb - grad_v.x* dt
		- grad_v.y* dt
		+ grad_v.z* dt;
	w_lattice[9] = w_nb - grad_w.x* dt
		- grad_w.y* dt
		+ grad_w.z* dt;




	///  case 11: back top west
	rho_lattice[11] = rho_nb - grad_rho.x* dt
		+ grad_rho.y* dt
		- grad_rho.z* dt;
	u_lattice[11] = u_nb - grad_u.x* dt
		+ grad_u.y* dt
		- grad_u.z* dt;
	v_lattice[11] = v_nb - grad_v.x* dt
		+ grad_v.y* dt
		- grad_v.z* dt;
	w_lattice[11] = w_nb - grad_w.x* dt
		+ grad_w.y* dt
		- grad_w.z* dt;




	/// case 14: front top west

	rho_lattice[14] = rho_nb - grad_rho.x* dt
		+ grad_rho.y* dt
		+ grad_rho.z* dt;
	u_lattice[14] = u_nb - grad_u.x* dt
		+ grad_u.y* dt
		+ grad_u.z* dt;
	v_lattice[14] = v_nb - grad_v.x* dt
		+ grad_v.y* dt
		+ grad_v.z* dt;
	w_lattice[14] = w_nb - grad_w.x* dt
		+ grad_w.y* dt
		+ grad_w.z* dt;




	/// case 8: front top east

	rho_lattice[8] = rho_nb + grad_rho.x* dt
		+ grad_rho.y* dt
		+ grad_rho.z* dt;
	u_lattice[8] = u_nb + grad_u.x* dt
		+ grad_u.y* dt
		+ grad_u.z* dt;
	v_lattice[8] = v_nb + grad_v.x* dt
		+ grad_v.y* dt
		+ grad_v.z* dt;
	w_lattice[8] = w_nb + grad_w.x* dt
		+ grad_w.y* dt
		+ grad_w.z* dt;



	/// case 10 Back Top East

	rho_lattice[10] = rho_nb + grad_rho.x* dt
		+ grad_rho.y* dt
		- grad_rho.z* dt;
	u_lattice[10] = u_nb + grad_u.x* dt
		+ grad_u.y* dt
		- grad_u.z* dt;
	v_lattice[10] = v_nb + grad_v.x* dt
		+ grad_v.y* dt
		- grad_v.z* dt;
	w_lattice[10] = w_nb + grad_w.x* dt
		+ grad_w.y* dt
		- grad_w.z* dt;


	/// case 12 Front Bottom East

	rho_lattice[12] = rho_nb + grad_rho.x* dt
		- grad_rho.y* dt
		+ grad_rho.z* dt;
	u_lattice[12] = u_nb + grad_u.x* dt
		- grad_u.y* dt
		+ grad_u.z* dt;
	v_lattice[12] = v_nb + grad_v.x* dt
		- grad_v.y* dt
		+ grad_v.z* dt;
	w_lattice[12] = w_nb + grad_w.x* dt
		- grad_w.y* dt
		+ grad_w.z* dt;


	/// case 13 Back Bottom East

	rho_lattice[13] = rho_nb + grad_rho.x* dt
		- grad_rho.y* dt
		- grad_rho.z* dt;
	u_lattice[13] = u_nb + grad_u.x* dt
		- grad_u.y* dt
		- grad_u.z* dt;
	v_lattice[13] = v_nb + grad_v.x* dt
		- grad_v.y* dt
		- grad_v.z* dt;
	w_lattice[13] = w_nb + grad_w.x* dt
		- grad_w.y* dt
		- grad_w.z* dt;


}



__device__ void populate_lattice_macros(double u_lattice[], double v_lattice[],
	double w_lattice[], double rho_lattice[], double3 cell_1, double3 cell_2,
	double3 interface_node, int i, int neighbour,
	double3 grad_u_1, double3 grad_u_2, double3 grad_v_1, double3  grad_v_2, double3  grad_w_1, double3  grad_w_2, double3  grad_rho_1, double3  grad_rho_2,
	double4 owner_soln, double4 neighbour_soln, int n_cells, double3 cell_normal, double dt, int bcs_rho_type[], int bcs_vel_type[], double4 bcs_macros[],
	double d1, double d173) {

	double3 temp1, temp2, temp3;

	///   case 0: // center node

	temp1 = cell_1 - interface_node;
	temp2 = cell_2 - interface_node;
	temp3 = cell_2 - cell_1;

	double rho_i, rho_nb, u_i, u_nb, v_i, v_nb, w_i, w_nb;
	double3 grad_rho_3, grad_u_3, grad_v_3, grad_w_3; // interpolated gradients for cell interface
	double c1, c2, c_mag, c_mag1, c_mag2;

	int nb;
	nb = neighbour - n_cells;
	/*double d_1 = d1 * dt;
	double d_173 = d173 * dt;*/
	double d_1 = 0.001 * dt;
	double d_173 = 0.001 * dt;

	if (neighbour > n_cells) {
		double4 bcs = bcs_macros[nb];
		// vboundary = v_i + grad_boundary * distance _ib

		if (bcs_rho_type[nb] == 1) {
			rho_lattice[0] = bcs.w;


		}
		else {
			rho_lattice[0] = owner_soln.w - dot_product(temp1, grad_rho_2);
		}

		if (bcs_vel_type[nb] == 1  ) {
			u_lattice[0] = bcs.x;
			v_lattice[0] = bcs.y;
			w_lattice[0] = bcs.z;


		}
		else {
			u_lattice[0] = owner_soln.x - dot_product(temp1, grad_u_2);
			v_lattice[0] = owner_soln.y - dot_product(temp1, grad_v_2);
			w_lattice[0] = owner_soln.z - dot_product(temp1, grad_w_2);
		}

		rho_i = rho_lattice[0];
		rho_nb = rho_lattice[0];

		u_i = u_lattice[0];
		u_nb = u_lattice[0];

		v_i = v_lattice[0];
		v_nb = v_lattice[0];

		w_i = w_lattice[0];
		w_nb = w_lattice[0];
		grad_rho_3 = grad_rho_2;
		grad_u_3 = grad_u_2;
		grad_v_3 = grad_v_2;
		grad_w_3 = grad_w_2;

		grad_rho_1 = grad_rho_2;
		grad_u_1 = grad_u_2;
		grad_v_1 = grad_v_2;
		grad_w_1 = grad_w_2;


	}
	else {

		rho_i = owner_soln.w - dot_product(temp1, grad_rho_1);
		rho_nb = neighbour_soln.w - dot_product(temp2, grad_rho_2);
		rho_lattice[0] = (rho_i + rho_nb)*0.5;

		u_i = owner_soln.x - dot_product(temp1, grad_u_1);
		u_nb = neighbour_soln.x - dot_product(temp2, grad_u_2);
		u_lattice[0] = (u_i + u_nb)*0.5;

		v_i = owner_soln.y - dot_product(temp1, grad_v_1);
		v_nb = neighbour_soln.y - dot_product(temp2, grad_v_2);
		v_lattice[0] = (v_i + v_nb)*0.5;

		w_i = owner_soln.z - dot_product(temp1, grad_w_1);
		w_nb = neighbour_soln.z - dot_product(temp2, grad_w_2);
		w_lattice[0] = (w_i + w_nb)*0.5;

		//interpolate gradients for cell interface gradients
		grad_u_3.x = (grad_u_1.x + grad_u_2.x)*0.5;
		grad_u_3.y = (grad_u_1.y + grad_u_2.y)*0.5;
		grad_u_3.z = (grad_u_1.z + grad_u_2.z)*0.5;

		grad_v_3.x = (grad_v_1.x + grad_v_2.x)*0.5;
		grad_v_3.y = (grad_v_1.y + grad_v_2.y)*0.5;
		grad_v_3.z = (grad_v_1.z + grad_v_2.z)*0.5;

		grad_w_3.x = (grad_w_1.x + grad_w_2.x)*0.5;
		grad_w_3.y = (grad_w_1.y + grad_w_2.y)*0.5;
		grad_w_3.z = (grad_w_1.z + grad_w_2.z)*0.5;

		grad_rho_3.x = (grad_rho_1.x + grad_rho_2.x)*0.5;
		grad_rho_3.y = (grad_rho_1.y + grad_rho_2.y)*0.5;
		grad_rho_3.z = (grad_rho_1.z + grad_rho_2.z)*0.5;

		/*grad_rho_3.x = c2 * grad_rho_1.x + c1 * grad_rho_2.x;
		grad_u_3.x = c2 * grad_rho_1.x + c1 * grad_rho_2.x;
		grad_v_3.x = c2 * grad_rho_1.x + c1 * grad_rho_2.x;
		grad_w_3.x = c2 * grad_rho_1.x + c1 * grad_rho_2.x;

		grad_rho_3.y = c2 * grad_rho_1.y + c1 * grad_rho_2.y;
		grad_u_3.y = c2 * grad_rho_1.y + c1 * grad_rho_2.y;
		grad_v_3.y = c2 * grad_rho_1.y + c1 * grad_rho_2.y;
		grad_w_3.y = c2 * grad_rho_1.y + c1 * grad_rho_2.y;

		grad_rho_3.z = c2 * grad_rho_1.z + c1 * grad_rho_2.z;
		grad_u_3.z = c2 * grad_rho_1.z + c1 * grad_rho_2.z;
		grad_v_3.z = c2 * grad_rho_1.z + c1 * grad_rho_2.z;
		grad_w_3.z = c2 * grad_rho_1.z + c1 * grad_rho_2.z;*/

	}


	///  case 1 and 2 :west_node and east_node

	if (cell_normal.x > d_1) {
		rho_lattice[1] = rho_i - grad_rho_1.x* dt;
		u_lattice[1] = u_i - grad_u_1.x* dt;
		v_lattice[1] = v_i - grad_v_1.x *dt;
		w_lattice[1] = w_i - grad_w_1.x *dt;

		rho_lattice[2] = rho_nb + grad_rho_2.x* dt;
		u_lattice[2] = u_nb + grad_u_2.x* dt;
		v_lattice[2] = v_nb + grad_v_2.x *dt;
		w_lattice[2] = w_nb + grad_w_2.x *dt;
	}
	else if (cell_normal.x < -d_1) {

		rho_lattice[1] = rho_nb - grad_rho_2.x* dt;
		u_lattice[1] = u_nb - grad_u_2.x* dt;
		v_lattice[1] = v_nb - grad_v_2.x *dt;
		w_lattice[1] = w_nb - grad_w_2.x *dt;

		rho_lattice[2] = rho_i + grad_rho_1.x* dt;
		u_lattice[2] = u_i + grad_u_1.x* dt;
		v_lattice[2] = v_i + grad_v_1.x *dt;
		w_lattice[2] = w_i + grad_w_1.x *dt;
	}
	else {
		rho_lattice[1] = rho_lattice[0] - grad_rho_3.x* dt;
		u_lattice[1] = u_lattice[0] - grad_u_3.x* dt;
		v_lattice[1] = v_lattice[0] - grad_v_3.x *dt;
		w_lattice[1] = w_lattice[0] - grad_w_3.x *dt;

		rho_lattice[2] = rho_lattice[0] + grad_rho_3.x* dt;
		u_lattice[2] = u_lattice[0] + grad_u_3.x* dt;
		v_lattice[2] = v_lattice[0] + grad_v_3.x *dt;
		w_lattice[2] = w_lattice[0] + grad_w_3.x *dt;

	}

	///   case 3 and 4: // bottom node and top node

	if (cell_normal.y > d_1) {

		rho_lattice[3] = rho_i - grad_rho_1.y* dt;
		u_lattice[3] = u_i - grad_u_1.y* dt;
		v_lattice[3] = v_i - grad_v_1.y *dt;
		w_lattice[3] = w_i - grad_w_1.y *dt;

		rho_lattice[4] = rho_nb + grad_rho_2.y* dt;
		u_lattice[4] = u_nb + grad_u_2.y* dt;
		v_lattice[4] = v_nb + grad_v_2.y *dt;
		w_lattice[4] = w_nb + grad_w_2.y *dt;


	}
	else if (cell_normal.y < -d_1) {
		rho_lattice[3] = rho_nb - grad_rho_2.y* dt;
		u_lattice[3] = u_nb - grad_u_2.y* dt;
		v_lattice[3] = v_nb - grad_v_2.y *dt;
		w_lattice[3] = w_nb - grad_w_2.y *dt;

		rho_lattice[4] = rho_i + grad_rho_1.y* dt;
		u_lattice[4] = u_i + grad_u_1.y* dt;
		v_lattice[4] = v_i + grad_v_1.y *dt;
		w_lattice[4] = w_i + grad_w_1.y *dt;

	}
	else {
		rho_lattice[3] = rho_lattice[0] - grad_rho_3.y* dt;
		u_lattice[3] = u_lattice[0] - grad_u_3.y* dt;
		v_lattice[3] = v_lattice[0] - grad_v_3.y *dt;
		w_lattice[3] = w_lattice[0] - grad_w_3.y *dt;

		rho_lattice[4] = rho_lattice[0] + grad_rho_3.y* dt;
		u_lattice[4] = u_lattice[0] + grad_u_3.y* dt;
		v_lattice[4] = v_lattice[0] + grad_v_3.y *dt;
		w_lattice[4] = w_lattice[0] + grad_w_3.y *dt;

	}




	///   case 5 and 6: // back node and front node

	if (cell_normal.z > d_1) {

		rho_lattice[5] = rho_i - grad_rho_1.z* dt;
		u_lattice[5] = u_i - grad_u_1.z* dt;
		v_lattice[5] = v_i - grad_v_1.z *dt;
		w_lattice[5] = w_i - grad_w_1.z *dt;

		rho_lattice[6] = rho_nb + grad_rho_2.z* dt;
		u_lattice[6] = u_nb + grad_u_2.z* dt;
		v_lattice[6] = v_nb + grad_v_2.z *dt;
		w_lattice[6] = w_nb + grad_w_2.z *dt;

	}
	else if (cell_normal.z < -d_1) {
		rho_lattice[5] = rho_nb - grad_rho_2.z* dt;
		u_lattice[5] = u_nb - grad_u_2.z* dt;
		v_lattice[5] = v_nb - grad_v_2.z *dt;
		w_lattice[5] = w_nb - grad_w_2.z *dt;

		rho_lattice[6] = rho_i + grad_rho_1.z* dt;
		u_lattice[6] = u_i + grad_u_1.z* dt;
		v_lattice[6] = v_i + grad_v_1.z *dt;
		w_lattice[6] = w_i + grad_w_1.z *dt;
	}
	else {
		rho_lattice[5] = rho_lattice[0] - grad_rho_3.z* dt;
		u_lattice[5] = u_lattice[0] - grad_u_3.z* dt;
		v_lattice[5] = v_lattice[0] - grad_v_3.z *dt;
		w_lattice[5] = w_lattice[0] - grad_w_3.z *dt;

		rho_lattice[6] = rho_lattice[0] + grad_rho_3.z* dt;
		u_lattice[6] = u_lattice[0] + grad_u_3.z* dt;
		v_lattice[6] = v_lattice[0] + grad_v_3.z *dt;
		w_lattice[6] = w_lattice[0] + grad_w_3.z *dt;

	}


	///  case 7: back bottom west and case 8: front top east
	if ((cell_normal.x + cell_normal.y + cell_normal.z) > d_173) {

		rho_lattice[7] = rho_i - grad_rho_1.x* dt
			- grad_rho_1.y* dt
			- grad_rho_1.z* dt;
		u_lattice[7] = u_i - grad_u_1.x* dt
			- grad_u_1.y* dt
			- grad_u_1.z* dt;
		v_lattice[7] = v_i - grad_v_1.x* dt
			- grad_v_1.y* dt
			- grad_v_1.z* dt;
		w_lattice[7] = w_i - grad_w_1.x* dt
			- grad_w_1.y* dt
			- grad_w_1.z* dt;

		rho_lattice[8] = rho_nb + grad_rho_2.x* dt
			+ grad_rho_2.y* dt
			+ grad_rho_2.z* dt;
		u_lattice[8] = u_nb + grad_u_2.x* dt
			+ grad_u_2.y* dt
			+ grad_u_2.z* dt;
		v_lattice[8] = v_nb + grad_v_2.x* dt
			+ grad_v_2.y* dt
			+ grad_v_2.z* dt;
		w_lattice[8] = w_nb + grad_w_2.x* dt
			+ grad_w_2.y* dt
			+ grad_w_2.z* dt;

	}
	else if ((cell_normal.x + cell_normal.y + cell_normal.z) < -d_173) {

		rho_lattice[7] = rho_nb - grad_rho_2.x* dt
			- grad_rho_2.y* dt
			- grad_rho_2.z* dt;
		u_lattice[7] = u_nb - grad_u_2.x* dt
			- grad_u_2.y* dt
			- grad_u_2.z* dt;
		v_lattice[7] = v_nb - grad_v_2.x* dt
			- grad_v_2.y* dt
			- grad_v_2.z* dt;
		w_lattice[7] = w_nb - grad_w_2.x* dt
			- grad_w_2.y* dt
			- grad_w_2.z* dt;


		rho_lattice[8] = rho_i + grad_rho_1.x* dt
			+ grad_rho_1.y* dt
			+ grad_rho_1.z* dt;
		u_lattice[8] = u_i + grad_u_1.x* dt
			+ grad_u_1.y* dt
			+ grad_u_1.z* dt;
		v_lattice[8] = v_i + grad_v_1.x* dt
			+ grad_v_1.y* dt
			+ grad_v_1.z* dt;
		w_lattice[8] = w_i + grad_w_1.x* dt
			+ grad_w_1.y* dt
			+ grad_w_1.z* dt;
	}
	else {
		rho_lattice[7] = rho_lattice[0] - grad_rho_3.x* dt
			- grad_rho_3.y* dt
			- grad_rho_3.z* dt;
		u_lattice[7] = u_lattice[0] - grad_u_3.x* dt
			- grad_u_3.y* dt
			- grad_u_3.z* dt;
		v_lattice[7] = v_lattice[0] - grad_v_3.x* dt
			- grad_v_3.y* dt
			- grad_v_3.z* dt;
		w_lattice[7] = w_lattice[0] - grad_w_3.x* dt
			- grad_w_3.y* dt
			- grad_w_3.z* dt;


		rho_lattice[8] = rho_lattice[0] + grad_rho_3.x* dt
			+ grad_rho_3.y* dt
			+ grad_rho_3.z* dt;
		u_lattice[8] = u_lattice[0] + grad_u_3.x* dt
			+ grad_u_3.y* dt
			+ grad_u_3.z* dt;
		v_lattice[8] = v_lattice[0] + grad_v_3.x* dt
			+ grad_v_3.y* dt
			+ grad_v_3.z* dt;
		w_lattice[8] = w_lattice[0] + grad_w_3.x* dt
			+ grad_w_3.y* dt
			+ grad_w_3.z* dt;

	}


	/// case 9: front bottom west and case 10 Back Top East
	if ((cell_normal.x + cell_normal.y + -1 * cell_normal.z) > d_173) {

		rho_lattice[9] = rho_i - grad_rho_1.x* dt
			- grad_rho_1.y* dt
			+ grad_rho_1.z* dt;
		u_lattice[9] = u_i - grad_u_1.x* dt
			- grad_u_1.y* dt
			+ grad_u_1.z* dt;
		v_lattice[9] = v_i - grad_v_1.x* dt
			- grad_v_1.y* dt
			+ grad_v_1.z* dt;
		w_lattice[9] = w_i - grad_w_1.x* dt
			- grad_w_1.y* dt
			+ grad_w_1.z* dt;

		rho_lattice[10] = rho_nb + grad_rho_2.x* dt
			+ grad_rho_2.y* dt
			- grad_rho_2.z* dt;
		u_lattice[10] = u_nb + grad_u_2.x* dt
			+ grad_u_2.y* dt
			- grad_u_2.z* dt;
		v_lattice[10] = v_nb + grad_v_2.x* dt
			+ grad_v_2.y* dt
			- grad_v_2.z* dt;
		w_lattice[10] = w_nb + grad_w_2.x* dt
			+ grad_w_2.y* dt
			- grad_w_2.z* dt;

	}
	else if ((cell_normal.x + cell_normal.y + -1 * cell_normal.z) < -d_173) {

		rho_lattice[9] = rho_nb - grad_rho_2.x* dt
			- grad_rho_2.y* dt
			+ grad_rho_2.z* dt;
		u_lattice[9] = u_nb - grad_u_2.x* dt
			- grad_u_2.y* dt
			+ grad_u_2.z* dt;
		v_lattice[9] = v_nb - grad_v_2.x* dt
			- grad_v_2.y* dt
			+ grad_v_2.z* dt;
		w_lattice[9] = w_nb - grad_w_2.x* dt
			- grad_w_2.y* dt
			+ grad_w_2.z* dt;

		rho_lattice[10] = rho_i + grad_rho_1.x* dt
			+ grad_rho_1.y* dt
			- grad_rho_1.z* dt;
		u_lattice[10] = u_i + grad_u_1.x* dt
			+ grad_u_1.y* dt
			- grad_u_1.z* dt;
		v_lattice[10] = v_i + grad_v_1.x* dt
			+ grad_v_1.y* dt
			- grad_v_1.z* dt;
		w_lattice[10] = w_i + grad_w_1.x* dt
			+ grad_w_1.y* dt
			- grad_w_1.z* dt;
	}
	else {

		rho_lattice[9] = rho_lattice[0] - grad_rho_3.x* dt
			- grad_rho_3.y* dt
			+ grad_rho_3.z* dt;
		u_lattice[9] = u_lattice[0] - grad_u_3.x* dt
			- grad_u_3.y* dt
			+ grad_u_3.z* dt;
		v_lattice[9] = v_lattice[0] - grad_v_3.x* dt
			- grad_v_3.y* dt
			+ grad_v_3.z* dt;
		w_lattice[9] = w_lattice[0] - grad_w_3.x* dt
			- grad_w_3.y* dt
			+ grad_w_3.z* dt;

		rho_lattice[10] = rho_lattice[0] + grad_rho_3.x* dt
			+ grad_rho_3.y* dt
			- grad_rho_3.z* dt;
		u_lattice[10] = u_lattice[0] + grad_u_3.x* dt
			+ grad_u_3.y* dt
			- grad_u_3.z* dt;
		v_lattice[10] = v_lattice[0] + grad_v_3.x* dt
			+ grad_v_3.y* dt
			- grad_v_3.z* dt;
		w_lattice[10] = w_lattice[0] + grad_w_3.x* dt
			+ grad_w_3.y* dt
			- grad_w_3.z* dt;


	}

	/// case 11: back top west and case 12 Front Bottom East
	if ((cell_normal.x + -1 * cell_normal.y + cell_normal.z) > d_173) {

		rho_lattice[11] = rho_i - grad_rho_1.x* dt
			+ grad_rho_1.y* dt
			- grad_rho_1.z* dt;
		u_lattice[11] = u_i - grad_u_1.x* dt
			+ grad_u_1.y* dt
			- grad_u_1.z* dt;
		v_lattice[11] = v_i - grad_v_1.x* dt
			+ grad_v_1.y* dt
			- grad_v_1.z* dt;
		w_lattice[11] = w_i - grad_w_1.x* dt
			+ grad_w_1.y* dt
			- grad_w_1.z* dt;

		rho_lattice[12] = rho_nb + grad_rho_2.x* dt
			- grad_rho_2.y* dt
			+ grad_rho_2.z* dt;
		u_lattice[12] = u_nb + grad_u_2.x* dt
			- grad_u_2.y* dt
			+ grad_u_2.z* dt;
		v_lattice[12] = v_nb + grad_v_2.x* dt
			- grad_v_2.y* dt
			+ grad_v_2.z* dt;
		w_lattice[12] = w_nb + grad_w_2.x* dt
			- grad_w_2.y* dt
			+ grad_w_2.z* dt;

	}
	else if ((cell_normal.x + -1 * cell_normal.y + cell_normal.z) < -d_173) {

		rho_lattice[11] = rho_nb - grad_rho_2.x* dt
			+ grad_rho_2.y* dt
			- grad_rho_2.z* dt;
		u_lattice[11] = u_nb - grad_u_2.x* dt
			+ grad_u_2.y* dt
			- grad_u_2.z* dt;
		v_lattice[11] = v_nb - grad_v_2.x* dt
			+ grad_v_2.y* dt
			- grad_v_2.z* dt;
		w_lattice[11] = w_nb - grad_w_2.x* dt
			+ grad_w_2.y* dt
			- grad_w_2.z* dt;

		rho_lattice[12] = rho_i + grad_rho_1.x* dt
			- grad_rho_1.y* dt
			+ grad_rho_1.z* dt;
		u_lattice[12] = u_i + grad_u_1.x* dt
			- grad_u_1.y* dt
			+ grad_u_1.z* dt;
		v_lattice[12] = v_i + grad_v_1.x* dt
			- grad_v_1.y* dt
			+ grad_v_1.z* dt;
		w_lattice[12] = w_i + grad_w_1.x* dt
			- grad_w_1.y* dt
			+ grad_w_1.z* dt;
	}
	else {
		rho_lattice[11] = rho_lattice[0] - grad_rho_3.x* dt
			+ grad_rho_3.y* dt
			- grad_rho_3.z* dt;
		u_lattice[11] = u_lattice[0] - grad_u_3.x* dt
			+ grad_u_3.y* dt
			- grad_u_3.z* dt;
		v_lattice[11] = v_lattice[0] - grad_v_3.x* dt
			+ grad_v_3.y* dt
			- grad_v_3.z* dt;
		w_lattice[11] = w_lattice[0] - grad_w_3.x* dt
			+ grad_w_3.y* dt
			- grad_w_3.z* dt;

		rho_lattice[12] = rho_lattice[0] + grad_rho_3.x* dt
			- grad_rho_3.y* dt
			+ grad_rho_3.z* dt;
		u_lattice[12] = u_lattice[0] + grad_u_3.x* dt
			- grad_u_3.y* dt
			+ grad_u_3.z* dt;
		v_lattice[12] = v_lattice[0] + grad_v_3.x* dt
			- grad_v_3.y* dt
			+ grad_v_3.z* dt;
		w_lattice[12] = w_lattice[0] + grad_w_3.x* dt
			- grad_w_3.y* dt
			+ grad_w_3.z* dt;
	}



	/// case 13 Back Bottom East and  case 14: front top west
	if ((-1 * cell_normal.x + cell_normal.y + cell_normal.z) > d_173) {

		rho_lattice[13] = rho_i + grad_rho_1.x* dt
			- grad_rho_1.y* dt
			- grad_rho_1.z* dt;
		u_lattice[13] = u_i + grad_u_1.x* dt
			- grad_u_1.y* dt
			- grad_u_1.z* dt;
		v_lattice[13] = v_i + grad_v_1.x* dt
			- grad_v_1.y* dt
			- grad_v_1.z* dt;
		w_lattice[13] = w_i + grad_w_1.x* dt
			- grad_w_1.y* dt
			- grad_w_1.z* dt;

		rho_lattice[14] = rho_nb - grad_rho_2.x* dt
			+ grad_rho_2.y* dt
			+ grad_rho_2.z* dt;
		u_lattice[14] = u_nb - grad_u_2.x* dt
			+ grad_u_2.y* dt
			+ grad_u_2.z* dt;
		v_lattice[14] = v_nb - grad_v_2.x* dt
			+ grad_v_2.y* dt
			+ grad_v_2.z* dt;
		w_lattice[14] = w_nb - grad_w_2.x* dt
			+ grad_w_2.y* dt
			+ grad_w_2.z* dt;

	}
	else if ((-1 * cell_normal.x + cell_normal.y + cell_normal.z) < -d_173) {

		rho_lattice[13] = rho_nb + grad_rho_2.x* dt
			- grad_rho_2.y* dt
			- grad_rho_2.z* dt;
		u_lattice[13] = u_nb + grad_u_2.x* dt
			- grad_u_2.y* dt
			- grad_u_2.z* dt;
		v_lattice[13] = v_nb + grad_v_2.x* dt
			- grad_v_2.y* dt
			- grad_v_2.z* dt;
		w_lattice[13] = w_nb + grad_w_2.x* dt
			- grad_w_2.y* dt
			- grad_w_2.z* dt;


		rho_lattice[14] = rho_i - grad_rho_1.x* dt
			+ grad_rho_1.y* dt
			+ grad_rho_1.z* dt;
		u_lattice[14] = u_i - grad_u_1.x* dt
			+ grad_u_1.y* dt
			+ grad_u_1.z* dt;
		v_lattice[14] = v_i - grad_v_1.x* dt
			+ grad_v_1.y* dt
			+ grad_v_1.z* dt;
		w_lattice[14] = w_i - grad_w_1.x* dt
			+ grad_w_1.y* dt
			+ grad_w_1.z* dt;
	}
	else {
		rho_lattice[13] = rho_lattice[0] + grad_rho_3.x* dt
			- grad_rho_3.y* dt
			- grad_rho_3.z* dt;
		u_lattice[13] = u_lattice[0] + grad_u_3.x* dt
			- grad_u_3.y* dt
			- grad_u_3.z* dt;
		v_lattice[13] = v_lattice[0] + grad_v_3.x* dt
			- grad_v_3.y* dt
			- grad_v_3.z* dt;
		w_lattice[13] = w_lattice[0] + grad_w_3.x* dt
			- grad_w_3.y* dt
			- grad_w_3.z* dt;


		rho_lattice[14] = rho_lattice[0] - grad_rho_3.x* dt
			+ grad_rho_3.y* dt
			+ grad_rho_3.z* dt;
		u_lattice[14] = u_lattice[0] - grad_u_3.x* dt
			+ grad_u_3.y* dt
			+ grad_u_3.z* dt;
		v_lattice[14] = v_lattice[0] - grad_v_3.x* dt
			+ grad_v_3.y* dt
			+ grad_v_3.z* dt;
		w_lattice[14] = w_lattice[0] - grad_w_3.x* dt
			+ grad_w_3.y* dt
			+ grad_w_3.z* dt;

	}



}



__device__ void populate_feq(double u_lattice[], double v_lattice[],
	double w_lattice[], double rho_lattice[],
	double feq_lattice[], int k, double pre_conditioned_gamma) {

	///d3q15 velocity set

	double uu2, vv2, u2v2w2, uv, uu, vv, ww2, uw, vw, ww;

	uu2 = u_lattice[k] * u_lattice[k] * pre_conditioned_gamma;
	vv2 = v_lattice[k] * v_lattice[k] * pre_conditioned_gamma;
	ww2 = w_lattice[k] * w_lattice[k] * pre_conditioned_gamma;
	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	uv = u_lattice[k] * v_lattice[k] * 9.0 * pre_conditioned_gamma;
	uw = u_lattice[k] * w_lattice[k] * 9.0 * pre_conditioned_gamma;
	vw = v_lattice[k] * w_lattice[k] * 9.0 * pre_conditioned_gamma;

	uu = u_lattice[k];
	vv = v_lattice[k];
	ww = w_lattice[k];

	switch (k) {

	case 0:
		feq_lattice[k] = lattice_weight[k] * rho_lattice[k] * (1.0

			- u2v2w2);
		break;

	case 1:
		feq_lattice[k] = lattice_weight[k] * rho_lattice[k] *
			(1.0 + 3.0*uu + 4.5*uu2 - u2v2w2);
		break;
	case 2:
		feq_lattice[k] = lattice_weight[k] * rho_lattice[k] *
			(1.0 - 3.0*uu + 4.5*uu2 - u2v2w2);
		break;
	case 3:
		feq_lattice[k] = lattice_weight[k] * rho_lattice[k] *
			(1.0 + 3.0*vv + 4.5*vv2 - u2v2w2);
		break;
	case 4:
		feq_lattice[k] = lattice_weight[k] * rho_lattice[k] *
			(1.0 - 3.0*vv + 4.5*vv2 - u2v2w2);
		break;
	case 5:
		feq_lattice[k] = lattice_weight[k] * rho_lattice[k] *
			(1.0 + 3.0*ww + 4.5*ww2 - u2v2w2);
		break;
	case 6:
		feq_lattice[k] = lattice_weight[k] * rho_lattice[k] *
			(1.0 - 3.0*ww + 4.5*ww2 - u2v2w2);
		break;
	case 7:
		feq_lattice[k] = lattice_weight[k] * rho_lattice[k] *
			(1.0 + 3.0*uu + 3.0*vv + 3.0*ww + uv + uw + vw +
				4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
		break;
	case 8:
		feq_lattice[k] = lattice_weight[k] * rho_lattice[k] *
			(1.0 - 3.0*uu - 3.0*vv - 3.0*ww + uv + uw + vw +
				4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
		break;
	case 9:
		feq_lattice[k] = lattice_weight[k] * rho_lattice[k] *
			(1.0 + 3.0*uu + 3.0*vv - 3.0*ww + uv - uw - vw +
				4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
		break;
	case 10:
		feq_lattice[k] = lattice_weight[k] * rho_lattice[k] *
			(1.0 - 3.0*uu - 3.0*vv + 3.0*ww + uv - uw - vw +
				4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
		break;
	case 11:
		feq_lattice[k] = lattice_weight[k] * rho_lattice[k] *
			(1.0 + 3.0*uu - 3.0*vv + 3.0*ww - uv + uw - vw +
				4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
		break;
	case 12:
		feq_lattice[k] = lattice_weight[k] * rho_lattice[k] *
			(1.0 - 3.0*uu + 3.0*vv - 3.0*ww - uv + uw - vw +
				4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
		break;
	case 13:
		feq_lattice[k] = lattice_weight[k] * rho_lattice[k] *
			(1.0 - 3.0*uu + 3.0*vv + 3.0*ww - uv - uw + vw +
				4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
		break;
	case 14:
		feq_lattice[k] = lattice_weight[k] * rho_lattice[k] *
			(1.0 + 3.0*uu - 3.0*vv - 3.0*ww - uv - uw + vw +
				4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
		break;
	}


}




__device__ void populate_feq_lattice_explicit(double u_lattice[], double v_lattice[],
	double w_lattice[], double rho_lattice[], double feq_lattice[]	) {

	///d3q15 velocity set
	int k;
	double uu2, vv2, u2v2w2, uv, uu, vv, ww2, uw, vw, ww;

	//trial with this value



	///////////////////////////////////////
	///K=0
	///////////////////////////////////////

	k = 0;
	uu2 = u_lattice[k] * u_lattice[k];
	vv2 = v_lattice[k] * v_lattice[k];
	ww2 = w_lattice[k] * w_lattice[k];
	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	feq_lattice[k ] = lattice_weight[k] * rho_lattice[k] * (1.0

		- u2v2w2);

	///////////////////////////////////////
	///K=1
	///////////////////////////////////////

	k = 1;
	uu2 = u_lattice[k] * u_lattice[k];
	vv2 = v_lattice[k] * v_lattice[k];
	ww2 = w_lattice[k] * w_lattice[k];
	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	

	uu = u_lattice[k];
	

	feq_lattice[k] = lattice_weight[k] * rho_lattice[k] *
		(1.0 + 3.0*uu + 4.5*uu2 - u2v2w2);

	///////////////////////////////////////
	///K=2
	///////////////////////////////////////

	k = 2;
	uu2 = u_lattice[k] * u_lattice[k];
	vv2 = v_lattice[k] * v_lattice[k];
	ww2 = w_lattice[k] * w_lattice[k];
	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	

	uu = u_lattice[k];
	

	feq_lattice[k] = lattice_weight[k] * rho_lattice[k] *
		(1.0 - 3.0*uu + 4.5*uu2 - u2v2w2);

	///////////////////////////////////////
	///K=3
	///////////////////////////////////////

	k = 3;
	uu2 = u_lattice[k] * u_lattice[k];
	vv2 = v_lattice[k] * v_lattice[k];
	ww2 = w_lattice[k] * w_lattice[k];
	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	

	
	vv = v_lattice[k];
	

	feq_lattice[ k] = lattice_weight[k] * rho_lattice[k] *
		(1.0 + 3.0*vv + 4.5*vv2 - u2v2w2);

	///////////////////////////////////////
	///K=4
	///////////////////////////////////////

	k = 4;
	uu2 = u_lattice[k] * u_lattice[k];
	vv2 = v_lattice[k] * v_lattice[k];
	ww2 = w_lattice[k] * w_lattice[k];
	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	

	vv = v_lattice[k];
	

	feq_lattice[ k] = lattice_weight[k] * rho_lattice[k] *
		(1.0 - 3.0*vv + 4.5*vv2 - u2v2w2);

	///////////////////////////////////////
	///K=5
	///////////////////////////////////////

	k = 5;
	uu2 = u_lattice[k] * u_lattice[k];
	vv2 = v_lattice[k] * v_lattice[k];
	ww2 = w_lattice[k] * w_lattice[k];
	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	

	
	ww = w_lattice[k];

	feq_lattice[ k] = lattice_weight[k] * rho_lattice[k] *
		(1.0 + 3.0*ww + 4.5*ww2 - u2v2w2);

	///////////////////////////////////////
	///K=6
	///////////////////////////////////////

	k = 6;
	uu2 = u_lattice[k] * u_lattice[k];
	vv2 = v_lattice[k] * v_lattice[k];
	ww2 = w_lattice[k] * w_lattice[k];
	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	

	
	ww = w_lattice[k];

	feq_lattice[ k] = lattice_weight[k] * rho_lattice[k] *
		(1.0 - 3.0*ww + 4.5*ww2 - u2v2w2);

	///////////////////////////////////////
	///K=7
	///////////////////////////////////////

	k = 7;
	uu2 = u_lattice[k] * u_lattice[k];
	vv2 = v_lattice[k] * v_lattice[k];
	ww2 = w_lattice[k] * w_lattice[k];
	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	uv = u_lattice[k] * v_lattice[k] * 9.0;
	uw = u_lattice[k] * w_lattice[k] * 9.0;
	vw = v_lattice[k] * w_lattice[k] * 9.0;

	uu = u_lattice[k];
	vv = v_lattice[k];
	ww = w_lattice[k];

	feq_lattice[ k] = lattice_weight[k] * rho_lattice[k] *
		(1.0 + 3.0*uu + 3.0*vv + 3.0*ww + uv + uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);

	///////////////////////////////////////
	///K=8
	///////////////////////////////////////

	k = 8;
	uu2 = u_lattice[k] * u_lattice[k];
	vv2 = v_lattice[k] * v_lattice[k];
	ww2 = w_lattice[k] * w_lattice[k];
	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	uv = u_lattice[k] * v_lattice[k] * 9.0;
	uw = u_lattice[k] * w_lattice[k] * 9.0;
	vw = v_lattice[k] * w_lattice[k] * 9.0;

	uu = u_lattice[k];
	vv = v_lattice[k];
	ww = w_lattice[k];

	feq_lattice[ k] = lattice_weight[k] * rho_lattice[k] *
		(1.0 - 3.0*uu - 3.0*vv - 3.0*ww + uv + uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);

	///////////////////////////////////////
	///K=9
	///////////////////////////////////////

	k = 9;
	uu2 = u_lattice[k] * u_lattice[k];
	vv2 = v_lattice[k] * v_lattice[k];
	ww2 = w_lattice[k] * w_lattice[k];
	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	uv = u_lattice[k] * v_lattice[k] * 9.0;
	uw = u_lattice[k] * w_lattice[k] * 9.0;
	vw = v_lattice[k] * w_lattice[k] * 9.0;

	uu = u_lattice[k];
	vv = v_lattice[k];
	ww = w_lattice[k];

	feq_lattice[ k] = lattice_weight[k] * rho_lattice[k] *
		(1.0 + 3.0*uu + 3.0*vv - 3.0*ww + uv - uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	///////////////////////////////////////
	///K=10
	///////////////////////////////////////

	k = 10;
	uu2 = u_lattice[k] * u_lattice[k];
	vv2 = v_lattice[k] * v_lattice[k];
	ww2 = w_lattice[k] * w_lattice[k];
	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	uv = u_lattice[k] * v_lattice[k] * 9.0;
	uw = u_lattice[k] * w_lattice[k] * 9.0;
	vw = v_lattice[k] * w_lattice[k] * 9.0;

	uu = u_lattice[k];
	vv = v_lattice[k];
	ww = w_lattice[k];

	feq_lattice[ k] = lattice_weight[k] * rho_lattice[k] *
		(1.0 - 3.0*uu - 3.0*vv + 3.0*ww + uv - uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);

	///////////////////////////////////////
	///K=11
	///////////////////////////////////////

	k = 11;
	uu2 = u_lattice[k] * u_lattice[k];
	vv2 = v_lattice[k] * v_lattice[k];
	ww2 = w_lattice[k] * w_lattice[k];
	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	uv = u_lattice[k] * v_lattice[k] * 9.0;
	uw = u_lattice[k] * w_lattice[k] * 9.0;
	vw = v_lattice[k] * w_lattice[k] * 9.0;

	uu = u_lattice[k];
	vv = v_lattice[k];
	ww = w_lattice[k];

	feq_lattice[ k] = lattice_weight[k] * rho_lattice[k] *
		(1.0 + 3.0*uu - 3.0*vv + 3.0*ww - uv + uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);

	///////////////////////////////////////
	///K=12
	///////////////////////////////////////

	k = 12;
	uu2 = u_lattice[k] * u_lattice[k];
	vv2 = v_lattice[k] * v_lattice[k];
	ww2 = w_lattice[k] * w_lattice[k];
	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	uv = u_lattice[k] * v_lattice[k] * 9.0;
	uw = u_lattice[k] * w_lattice[k] * 9.0;
	vw = v_lattice[k] * w_lattice[k] * 9.0;

	uu = u_lattice[k];
	vv = v_lattice[k];
	ww = w_lattice[k];

	feq_lattice[ k] = lattice_weight[k] * rho_lattice[k] *
		(1.0 - 3.0*uu + 3.0*vv - 3.0*ww - uv + uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);

	///////////////////////////////////////
	///K=13
	///////////////////////////////////////

	k = 13;
	uu2 = u_lattice[k] * u_lattice[k];
	vv2 = v_lattice[k] * v_lattice[k];
	ww2 = w_lattice[k] * w_lattice[k];
	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	uv = u_lattice[k] * v_lattice[k] * 9.0;
	uw = u_lattice[k] * w_lattice[k] * 9.0;
	vw = v_lattice[k] * w_lattice[k] * 9.0;

	uu = u_lattice[k];
	vv = v_lattice[k];
	ww = w_lattice[k];

	feq_lattice[ k] = lattice_weight[k] * rho_lattice[k] *
		(1.0 - 3.0*uu + 3.0*vv + 3.0*ww - uv - uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);

	///////////////////////////////////////
	///K=14
	///////////////////////////////////////

	k = 14;
	uu2 = u_lattice[k] * u_lattice[k];
	vv2 = v_lattice[k] * v_lattice[k];
	ww2 = w_lattice[k] * w_lattice[k];
	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	uv = u_lattice[k] * v_lattice[k] * 9.0;
	uw = u_lattice[k] * w_lattice[k] * 9.0;
	vw = v_lattice[k] * w_lattice[k] * 9.0;

	uu = u_lattice[k];
	vv = v_lattice[k];
	ww = w_lattice[k];

	feq_lattice[ k] = lattice_weight[k] * rho_lattice[k] *
		(1.0 + 3.0*uu - 3.0*vv - 3.0*ww - uv - uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);


	


}



__device__ void populate_feq_lattice(double u_lattice[], double v_lattice[],
	double w_lattice[], double rho_lattice[], double * feq_lattice,
	 int k,   int face) {

	///d3q15 velocity set

	double uu2, vv2, u2v2w2, uv, uu, vv, ww2, uw, vw, ww;

	//trial with this value

	int face_index = face * 15;


	uu2 = u_lattice[k] * u_lattice[k] ;
	vv2 = v_lattice[k] * v_lattice[k]  ;
	ww2 = w_lattice[k] * w_lattice[k]  ;
	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	uv = u_lattice[k] * v_lattice[k] * 9.0  ;
	uw = u_lattice[k] * w_lattice[k] * 9.0  ;
	vw = v_lattice[k] * w_lattice[k] * 9.0  ;

	uu = u_lattice[k];
	vv = v_lattice[k];
	ww = w_lattice[k];

	switch (k) {

	case 0:
		feq_lattice[face_index] = lattice_weight[k] * rho_lattice[k] * (1.0

			- u2v2w2);
		break;

	case 1:
		feq_lattice[face_index] = lattice_weight[k] * rho_lattice[k] *
			(1.0 + 3.0*uu + 4.5*uu2 - u2v2w2);
		break;
	case 2:
		feq_lattice[face_index] = lattice_weight[k] * rho_lattice[k] *
			(1.0 - 3.0*uu + 4.5*uu2 - u2v2w2);
		break;
	case 3:
		feq_lattice[face_index] = lattice_weight[k] * rho_lattice[k] *
			(1.0 + 3.0*vv + 4.5*vv2 - u2v2w2);
		break;
	case 4:
		feq_lattice[face_index] = lattice_weight[k] * rho_lattice[k] *
			(1.0 - 3.0*vv + 4.5*vv2 - u2v2w2);
		break;
	case 5:
		feq_lattice[face_index] = lattice_weight[k] * rho_lattice[k] *
			(1.0 + 3.0*ww + 4.5*ww2 - u2v2w2);
		break;
	case 6:
		feq_lattice[face_index] = lattice_weight[k] * rho_lattice[k] *
			(1.0 - 3.0*ww + 4.5*ww2 - u2v2w2);
		break;
	case 7:
		feq_lattice[face_index] = lattice_weight[k] * rho_lattice[k] *
			(1.0 + 3.0*uu + 3.0*vv + 3.0*ww + uv + uw + vw +
				4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
		break;
	case 8:
		feq_lattice[face_index] = lattice_weight[k] * rho_lattice[k] *
			(1.0 - 3.0*uu - 3.0*vv - 3.0*ww + uv + uw + vw +
				4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
		break;
	case 9:
		feq_lattice[face_index] = lattice_weight[k] * rho_lattice[k] *
			(1.0 + 3.0*uu + 3.0*vv - 3.0*ww + uv - uw - vw +
				4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
		break;
	case 10:
		feq_lattice[face_index] = lattice_weight[k] * rho_lattice[k] *
			(1.0 - 3.0*uu - 3.0*vv + 3.0*ww + uv - uw - vw +
				4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
		break;
	case 11:
		feq_lattice[face_index] = lattice_weight[k] * rho_lattice[k] *
			(1.0 + 3.0*uu - 3.0*vv + 3.0*ww - uv + uw - vw +
				4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
		break;
	case 12:
		feq_lattice[face_index] = lattice_weight[k] * rho_lattice[k] *
			(1.0 - 3.0*uu + 3.0*vv - 3.0*ww - uv + uw - vw +
				4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
		break;
	case 13:
		feq_lattice[face_index] = lattice_weight[k] * rho_lattice[k] *
			(1.0 - 3.0*uu + 3.0*vv + 3.0*ww - uv - uw + vw +
				4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
		break;
	case 14:
		feq_lattice[face_index] = lattice_weight[k] * rho_lattice[k] *
			(1.0 + 3.0*uu - 3.0*vv - 3.0*ww - uv - uw + vw +
				4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
		break;
	}


}


__device__ void calculate_flux_at_interface(double3 u_interface, double dt, double pre_conditioned_gamma,
	double rho_interface, double4 &cell_flux, int i,
	double* feq_lattice, double3 cell_normal, double interface_area, double* local_fneq, double visc, int testcase, double  visc_factor) {

	double uu2, vv2, ww2, u2v2w2, uu, vv, ww, uv, uw, vw, fneq_tau;
	double feq_interface[15];
	double4 x_flux, y_flux, z_flux;

	uu2 = u_interface.x * u_interface.x * pre_conditioned_gamma;
	vv2 = u_interface.y * u_interface.y * pre_conditioned_gamma;
	ww2 = u_interface.z *u_interface.z * pre_conditioned_gamma;

	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	uu = u_interface.x;
	vv = u_interface.y;
	ww = u_interface.z;

	uv = uu * vv*9.0 * pre_conditioned_gamma;
	uw = uu * ww*9.0 * pre_conditioned_gamma;
	vw = vv * ww *9.0 * pre_conditioned_gamma;

	if (testcase == 8) {
		fneq_tau = 0.0;
	}
	else {
		fneq_tau = (visc* visc_factor * 3 / dt * pre_conditioned_gamma);
	}
	
	local_fneq[i] = fneq_tau;


	feq_interface[1] = lattice_weight[1] * rho_interface*
		(1.0 + 3.0*uu + 4.5*uu2 - u2v2w2);
	feq_interface[1] = feq_interface[1]
		- fneq_tau * (feq_interface[1] - feq_lattice[1]);

	feq_interface[2] = lattice_weight[2] * rho_interface*
		(1.0 - 3.0*uu + 4.5*uu2 - u2v2w2);
	feq_interface[2] = feq_interface[2]
		- fneq_tau * (feq_interface[2] - feq_lattice[2]);

	feq_interface[3] = lattice_weight[3] * rho_interface*
		(1.0 + 3.0*vv + 4.5*vv2 - u2v2w2);
	feq_interface[3] = feq_interface[3]
		- fneq_tau * (feq_interface[3] - feq_lattice[3]);

	feq_interface[4] = lattice_weight[4] * rho_interface*
		(1.0 - 3.0*vv + 4.5*vv2 - u2v2w2);
	feq_interface[4] = feq_interface[4]
		- fneq_tau * (feq_interface[4] - feq_lattice[4]);

	feq_interface[5] = lattice_weight[5] * rho_interface*
		(1.0 + 3.0*ww + 4.5*ww2 - u2v2w2);
	feq_interface[5] = feq_interface[5]
		- fneq_tau * (feq_interface[5] - feq_lattice[5]);

	feq_interface[6] = lattice_weight[6] * rho_interface*
		(1.0 - 3.0*ww + 4.5*ww2 - u2v2w2);
	feq_interface[6] = feq_interface[6]
		- fneq_tau * (feq_interface[6] - feq_lattice[6]);

	feq_interface[7] = lattice_weight[7] * rho_interface*
		(1.0 + 3.0*uu + 3.0*vv + 3.0*ww + uv + uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[7] = feq_interface[7]
		- fneq_tau * (feq_interface[7] - feq_lattice[7]);

	feq_interface[8] = lattice_weight[8] * rho_interface*
		(1.0 - 3.0*uu - 3.0*vv - 3.0*ww + uv + uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[8] = feq_interface[8]
		- fneq_tau * (feq_interface[8] - feq_lattice[8]);

	feq_interface[9] = lattice_weight[9] * rho_interface*
		(1.0 + 3.0*uu + 3.0*vv - 3.0*ww + uv - uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[9] = feq_interface[9]
		- fneq_tau * (feq_interface[9] - feq_lattice[9]);

	feq_interface[10] = lattice_weight[10] * rho_interface*
		(1.0 - 3.0*uu - 3.0*vv + 3.0*ww + uv - uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[10] = feq_interface[10]
		- fneq_tau * (feq_interface[10] - feq_lattice[10]);

	feq_interface[11] = lattice_weight[11] * rho_interface*
		(1.0 + 3.0*uu - 3.0*vv + 3.0*ww - uv + uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[11] = feq_interface[11]
		- fneq_tau * (feq_interface[11] - feq_lattice[11]);

	feq_interface[12] = lattice_weight[12] * rho_interface*
		(1.0 - 3.0*uu + 3.0*vv - 3.0*ww - uv + uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[12] = feq_interface[12]
		- fneq_tau * (feq_interface[12] - feq_lattice[12]);

	feq_interface[13] = lattice_weight[13] * rho_interface*
		(1.0 - 3.0*uu + 3.0*vv + 3.0*ww - uv - uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[13] = feq_interface[13]
		- fneq_tau * (feq_interface[13] - feq_lattice[13]);

	feq_interface[14] = lattice_weight[14] * rho_interface*
		(1.0 + 3.0*uu - 3.0*vv - 3.0*ww - uv - uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[14] = feq_interface[14]
		- fneq_tau * (feq_interface[14] - feq_lattice[14]);


	x_flux.w = (feq_interface[1] - feq_interface[2]
		+ feq_interface[7] - feq_interface[8]
		+ feq_interface[9] - feq_interface[10] + feq_interface[11]
		- feq_interface[12] - feq_interface[13]
		+ feq_interface[14]);

	x_flux.x =
		(feq_interface[1] + feq_interface[2]
			+ feq_interface[7] + feq_interface[8]
			+ feq_interface[9] + feq_interface[10] + feq_interface[11]
			+ feq_interface[12] + feq_interface[13]
			+ feq_interface[14]);

	x_flux.y =
		feq_interface[7] + feq_interface[8]
		+ feq_interface[9] + feq_interface[10] - feq_interface[11] - feq_interface[12]
		- feq_interface[13]
		- feq_interface[14];

	x_flux.z =
		feq_interface[7] + feq_interface[8]
		- feq_interface[9] - feq_interface[10] + feq_interface[11] + feq_interface[12]
		- feq_interface[13]
		- feq_interface[14];



	y_flux.w = (feq_interface[3] - feq_interface[4]
		+ feq_interface[7] - feq_interface[8]
		+ feq_interface[9] - feq_interface[10] - feq_interface[11]
		+ feq_interface[12] + feq_interface[13]
		- feq_interface[14]);


	y_flux.x = x_flux.y;

	y_flux.y = (feq_interface[3] + feq_interface[4]
		+ feq_interface[7] + feq_interface[8]
		+ feq_interface[9] + feq_interface[10] + feq_interface[11]
		+ feq_interface[12] + feq_interface[13]
		+ feq_interface[14]);

	y_flux.z = feq_interface[7] + feq_interface[8]
		- feq_interface[9] - feq_interface[10] - feq_interface[11] - feq_interface[12]
		+ feq_interface[13]
		+ feq_interface[14];

	z_flux.w = (feq_interface[5] - feq_interface[6]
		+ feq_interface[7] - feq_interface[8]
		- feq_interface[9] + feq_interface[10] + feq_interface[11]
		- feq_interface[12] + feq_interface[13]
		- feq_interface[14]);
	z_flux.x = x_flux.z;
	z_flux.y = y_flux.z;
	z_flux.z = (feq_interface[5] + feq_interface[6]
		+ feq_interface[7] + feq_interface[8]
		+ feq_interface[9] + feq_interface[10] + feq_interface[11]
		+ feq_interface[12] + feq_interface[13]
		+ feq_interface[14]);


	cell_flux.w = (x_flux.w *cell_normal.x + y_flux.w* cell_normal.y +
		z_flux.w *cell_normal.z)*interface_area;
	cell_flux.x = (x_flux.x *cell_normal.x +
		y_flux.x * cell_normal.y +
		z_flux.x*cell_normal.z)*interface_area;

	cell_flux.y = (x_flux.y *cell_normal.x +
		y_flux.y * cell_normal.y +
		z_flux.y*cell_normal.z)*interface_area;


	cell_flux.z = (x_flux.z *cell_normal.x +
		y_flux.z * cell_normal.y +
		z_flux.z*cell_normal.z)*interface_area;
}



__device__ void calculate_flux_at_interface_incompressible(double3 u_interface, double dt, double pre_conditioned_gamma,
	double rho_interface, double4 &cell_flux, int i,
	double* feq_lattice, double3 cell_normal, double interface_area, double* local_fneq, double visc, int testcase, double  visc_factor) {

	double uu2, vv2, ww2, u2v2w2, uu, vv, ww, uv, uw, vw, fneq_tau;
	double feq_interface[15];
	double4 x_flux, y_flux, z_flux;

	uu2 = u_interface.x * u_interface.x * pre_conditioned_gamma;
	vv2 = u_interface.y * u_interface.y * pre_conditioned_gamma;
	ww2 = u_interface.z *u_interface.z * pre_conditioned_gamma;

	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	uu = u_interface.x;
	vv = u_interface.y;
	ww = u_interface.z;

	uv = uu * vv*9.0 * pre_conditioned_gamma;
	uw = uu * ww*9.0 * pre_conditioned_gamma;
	vw = vv * ww *9.0 * pre_conditioned_gamma;

	if (testcase == 8) {
		fneq_tau = 0.0;
	}
	else {
		fneq_tau = (visc* visc_factor * 3 / dt * pre_conditioned_gamma);
	}

	local_fneq[i] = fneq_tau;


	feq_interface[1] = lattice_weight[1] * (rho_interface +
		( 3.0*uu + 4.5*uu2 - u2v2w2));
	feq_interface[1] = feq_interface[1]
		- fneq_tau * (feq_interface[1] - feq_lattice[1]);

	feq_interface[2] = lattice_weight[2] * (rho_interface +
		( 3.0*uu + 4.5*uu2 - u2v2w2));
	feq_interface[2] = feq_interface[2]
		- fneq_tau * (feq_interface[2] - feq_lattice[2]);

	feq_interface[3] = lattice_weight[3] * (rho_interface +
		( 3.0*vv + 4.5*vv2 - u2v2w2));
	feq_interface[3] = feq_interface[3]
		- fneq_tau * (feq_interface[3] - feq_lattice[3]);

	feq_interface[4] = lattice_weight[4] * (rho_interface +
		( 3.0*vv + 4.5*vv2 - u2v2w2));
	feq_interface[4] = feq_interface[4]
		- fneq_tau * (feq_interface[4] - feq_lattice[4]);

	feq_interface[5] = lattice_weight[5] * (rho_interface +
		( 3.0*ww + 4.5*ww2 - u2v2w2));
	feq_interface[5] = feq_interface[5]
		- fneq_tau * (feq_interface[5] - feq_lattice[5]);

	feq_interface[6] = lattice_weight[6] * (rho_interface +
		( 3.0*ww + 4.5*ww2 - u2v2w2));
	feq_interface[6] = feq_interface[6]
		- fneq_tau * (feq_interface[6] - feq_lattice[6]);

	feq_interface[7] = lattice_weight[7] * (rho_interface +
		( 3.0*uu + 3.0*vv + 3.0*ww + uv + uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2));
	feq_interface[7] = feq_interface[7]
		- fneq_tau * (feq_interface[7] - feq_lattice[7]);

	feq_interface[8] = lattice_weight[8] * (rho_interface +
		( 3.0*uu - 3.0*vv - 3.0*ww + uv + uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2));
	feq_interface[8] = feq_interface[8]
		- fneq_tau * (feq_interface[8] - feq_lattice[8]);

	feq_interface[9] = lattice_weight[9] * (rho_interface +
		( 3.0*uu + 3.0*vv - 3.0*ww + uv - uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2));
	feq_interface[9] = feq_interface[9]
		- fneq_tau * (feq_interface[9] - feq_lattice[9]);

	feq_interface[10] = lattice_weight[10] * (rho_interface +
		( 3.0*uu - 3.0*vv + 3.0*ww + uv - uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2));
	feq_interface[10] = feq_interface[10]
		- fneq_tau * (feq_interface[10] - feq_lattice[10]);

	feq_interface[11] = lattice_weight[11] * (rho_interface +
		( 3.0*uu - 3.0*vv + 3.0*ww - uv + uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2));
	feq_interface[11] = feq_interface[11]
		- fneq_tau * (feq_interface[11] - feq_lattice[11]);

	feq_interface[12] = lattice_weight[12] * (rho_interface +
		( 3.0*uu + 3.0*vv - 3.0*ww - uv + uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2));
	feq_interface[12] = feq_interface[12]
		- fneq_tau * (feq_interface[12] - feq_lattice[12]);

	feq_interface[13] = lattice_weight[13] * (rho_interface +
		( 3.0*uu + 3.0*vv + 3.0*ww - uv - uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2));
	feq_interface[13] = feq_interface[13]
		- fneq_tau * (feq_interface[13] - feq_lattice[13]);

	feq_interface[14] = lattice_weight[14] * (rho_interface +
		( 3.0*uu - 3.0*vv - 3.0*ww - uv - uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2));
	feq_interface[14] = feq_interface[14]
		- fneq_tau * (feq_interface[14] - feq_lattice[14]);


	x_flux.w = (feq_interface[1] - feq_interface[2]
		+ feq_interface[7] - feq_interface[8]
		+ feq_interface[9] - feq_interface[10] + feq_interface[11]
		- feq_interface[12] - feq_interface[13]
		+ feq_interface[14]);

	x_flux.x =
		(feq_interface[1] + feq_interface[2]
			+ feq_interface[7] + feq_interface[8]
			+ feq_interface[9] + feq_interface[10] + feq_interface[11]
			+ feq_interface[12] + feq_interface[13]
			+ feq_interface[14]);

	x_flux.y =
		feq_interface[7] + feq_interface[8]
		+ feq_interface[9] + feq_interface[10] - feq_interface[11] - feq_interface[12]
		- feq_interface[13]
		- feq_interface[14];

	x_flux.z =
		feq_interface[7] + feq_interface[8]
		- feq_interface[9] - feq_interface[10] + feq_interface[11] + feq_interface[12]
		- feq_interface[13]
		- feq_interface[14];



	y_flux.w = (feq_interface[3] - feq_interface[4]
		+ feq_interface[7] - feq_interface[8]
		+ feq_interface[9] - feq_interface[10] - feq_interface[11]
		+ feq_interface[12] + feq_interface[13]
		- feq_interface[14]);


	y_flux.x = x_flux.y;

	y_flux.y = (feq_interface[3] + feq_interface[4]
		+ feq_interface[7] + feq_interface[8]
		+ feq_interface[9] + feq_interface[10] + feq_interface[11]
		+ feq_interface[12] + feq_interface[13]
		+ feq_interface[14]);

	y_flux.z = feq_interface[7] + feq_interface[8]
		- feq_interface[9] - feq_interface[10] - feq_interface[11] - feq_interface[12]
		+ feq_interface[13]
		+ feq_interface[14];

	z_flux.w = (feq_interface[5] - feq_interface[6]
		+ feq_interface[7] - feq_interface[8]
		- feq_interface[9] + feq_interface[10] + feq_interface[11]
		- feq_interface[12] + feq_interface[13]
		- feq_interface[14]);
	z_flux.x = x_flux.z;
	z_flux.y = y_flux.z;
	z_flux.z = (feq_interface[5] + feq_interface[6]
		+ feq_interface[7] + feq_interface[8]
		+ feq_interface[9] + feq_interface[10] + feq_interface[11]
		+ feq_interface[12] + feq_interface[13]
		+ feq_interface[14]);


	cell_flux.w = (x_flux.w *cell_normal.x + y_flux.w* cell_normal.y +
		z_flux.w *cell_normal.z)*interface_area;
	cell_flux.x = (x_flux.x *cell_normal.x +
		y_flux.x * cell_normal.y +
		z_flux.x*cell_normal.z)*interface_area;

	cell_flux.y = (x_flux.y *cell_normal.x +
		y_flux.y * cell_normal.y +
		z_flux.y*cell_normal.z)*interface_area;


	cell_flux.z = (x_flux.z *cell_normal.x +
		y_flux.z * cell_normal.y +
		z_flux.z*cell_normal.z)*interface_area;
}

// update bc nodes to allow for changes in solution

/// bc boundary conditions only ported for dirichlet and neumann conditions
__global__ void update_unstructured_bcs(int n_bc_cells, int n_neighbours, int n_cells, int* mesh_owner, int* bcs_rho_type, int* bcs_vel_type, double4* input_soln,
	double4* bc_arr, double3 * cell_centroid, double channel_diameter) {



	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	__shared__ double H2, H4;
	H2 = channel_diameter * 0.5;
	H4 = channel_diameter * channel_diameter *channel_diameter*channel_diameter;

	for (int i = index; i < n_bc_cells; i += stride) {

		if (i < n_bc_cells) {

			int j, face, nb;

			j = i + n_cells;
			face = i + n_neighbours;
			nb = mesh_owner[face];

			double4 nb_soln = input_soln[nb];
			double4 temp;
			double4 bc = bc_arr[i];


			///NEEDS to be modified for non-uniform solver
			// 1 = dirichlet, 2 = neumann, 3 = periodic
			if (bcs_rho_type[i] == 1) {
				temp.w = bc.w - (nb_soln.w - bc.w);
			}
			else if (bcs_rho_type[i] == 2) {

				temp.w = nb_soln.w;
			}

			if (bcs_vel_type[i] == 1) {

				temp.x = bc.x - (nb_soln.x - bc.x);
				temp.y = bc.y - (nb_soln.y - bc.y);
				temp.z = bc.z - (nb_soln.z - bc.z);
			}
			else if (bcs_vel_type[i] == 2) {

				temp.x = nb_soln.x;
				temp.y = nb_soln.y;
				temp.z = nb_soln.z;
			}
			else if (bcs_vel_type[i] == 9) {
				double3 cell = cell_centroid[j];
				double u_target;


				/// 4u ( H/2 +y) *(H/2 -y)/H^2
			/*	u_target = 4 * bc.x*(channel_diameter*0.2 / 0.41 + cell.y) * (channel_diameter - channel_diameter * 0.2 / 0.41 - cell.y) / (channel_diameter *channel_diameter);*/

				/// 16 u_max ( H/2 +y) *(H/2 -y) ( H/2 +z) *(H/2 -z)/H^4
				u_target = 16 * bc.x * (H2 + cell.y) *(H2 - cell.y) * (H2 + cell.z)*(H2 - cell.z) / H4;

				temp.x = u_target - (nb_soln.x - u_target);
				temp.y = bc.y - (nb_soln.y - bc.y);
				temp.z = bc.z - (nb_soln.z - bc.z);
			}
			else if (bcs_vel_type[i] == 10) { /// shear flow inlet
				double3 cell = cell_centroid[j];
				double u_target;


			
				// u_target =  (H2 + cell.y)*bc.x;
				u_target = (cell.y)*bc.x;
				temp.x = u_target - (nb_soln.x - u_target);
				temp.y = bc.y - (nb_soln.y - bc.y);
				temp.z = bc.z - (nb_soln.z - bc.z);
			}

			input_soln[j] = temp;

		}

	}



}



/// bc boundary conditions only ported for dirichlet and neumann conditions
__global__ void get_bc_gradients(int n_bc_cells, int n_neighbours, int n_cells, int* mesh_owner, int* bcs_rho_type, int* bcs_vel_type, double4* input_soln,
	double3* face_normal_arr, double3* centroid, double4* bc_arr,
	double3* grad_rho_arr, double3* grad_u_arr, double3* grad_v_arr, double3* grad_w_arr) {



	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	for (int i = index; i < n_bc_cells; i += stride) {

		if (i < n_bc_cells) {

			double3 cell_to_face = make_double3(0.0, 0.0, 0.0);
			double3 temp, neighbour, face_normal;


			int face = n_neighbours + i;
			int neighbour_id = mesh_owner[face];
			//periodic nodes get LS treatment
			int j = i + n_cells;

			face_normal = face_normal_arr[face];

			double3 centroid_nb = centroid[neighbour_id];
			double3 centroid_bc = centroid[j];
			double dx, dy, dz, d_mag, d_u, d_v, d_w, d_rho;
			double ni, nj, nk;

			double4 bc = input_soln[j];
			double4 soln = input_soln[neighbour_id];

			dx = centroid_bc.x - centroid_nb.x;
			dy = centroid_bc.y - centroid_nb.y;
			dz = centroid_bc.z - centroid_nb.z;

			d_mag = sqrt(pow(dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0));

			//get unit vectors
			ni = dx / d_mag;
			nj = dy / d_mag;
			nk = dz / d_mag;

			//bc is now at ghost cell centre, change
		/*	d_u = 2 * (bc.x - soln.x);
			d_v = 2 * (bc.y - soln.y);
			d_w = 2 * (bc.z - soln.z);
			d_rho = 2 * (bc.w - soln.w);*/

			d_u = (bc.x - soln.x);
			d_v =  (bc.y - soln.y);
			d_w = (bc.z - soln.z);
			d_rho =  (bc.w - soln.w);

			//dirichlet
			if (bcs_rho_type[i] == 1) {

				temp.x = d_rho / d_mag * ni;
				temp.y = d_rho / d_mag * nj;
				temp.z = d_rho / d_mag * nk;

				//neumann
			}
			else if (bcs_rho_type[i] == 2) {
				neighbour = grad_rho_arr[neighbour_id];
				temp.x = neighbour.x *(1 - fabs(face_normal.x));
				temp.y = neighbour.y*(1 - fabs(face_normal.y));
				temp.z = neighbour.z*(1 - fabs(face_normal.z));

			}
			grad_rho_arr[j] = temp;


			//dirichlet
			if (bcs_vel_type[i] == 1 || bcs_vel_type[i] == 9 || bcs_vel_type[i] == 10) {
				temp.x = d_u / d_mag * ni;
				temp.y = d_u / d_mag * nj;
				temp.z = d_u / d_mag * nk;
				grad_u_arr[j] = temp;

				temp.x = d_v / d_mag * ni;
				temp.y = d_v / d_mag * nj;
				temp.z = d_v / d_mag * nk;
				grad_v_arr[j] = temp;

				temp.x = d_w / d_mag * ni;
				temp.y = d_w / d_mag * nj;
				temp.z = d_w / d_mag * nk;
				grad_w_arr[j] = temp;
			}
			else if (bcs_vel_type[i] == 2) {

				// page 612 of Moukallad
				neighbour = grad_u_arr[neighbour_id];
				temp.x = neighbour.x *(1 - fabs(face_normal.x));
				temp.y = neighbour.y*(1 - fabs(face_normal.y));
				temp.z = neighbour.z*(1 - fabs(face_normal.z));
				grad_u_arr[j] = temp;

				neighbour = grad_v_arr[neighbour_id];
				temp.x = neighbour.x *(1 - fabs(face_normal.x));
				temp.y = neighbour.y*(1 - fabs(face_normal.y));
				temp.z = neighbour.z*(1 - fabs(face_normal.z));
				grad_v_arr[j] = temp;

				neighbour = grad_w_arr[neighbour_id];
				temp.x = neighbour.x *(1 - fabs(face_normal.x));
				temp.y = neighbour.y*(1 - fabs(face_normal.y));
				temp.z = neighbour.z*(1 - fabs(face_normal.z));
				grad_w_arr[j] = temp;

			}

		}
	}
}


__device__ void add_LS_contributions(double3 &RHS_rho, double3 &RHS_u, double3 &RHS_v, double3 &RHS_w, double4* src, int i1, int i, int i_nb, double3* RHS_array) {

	double d_u, d_v, d_w, d_rho;

	//boundary condition, find gradient at shared face
	double4 cell_1, cell_2;
	cell_1 = src[i1];
	cell_2 = src[i];
	double3 RHS;
	RHS = RHS_array[i * 6 + i_nb];

	d_u = (cell_1.x - cell_2.x);
	d_v = (cell_1.y - cell_2.y);
	d_w = (cell_1.z - cell_2.z);
	d_rho = (cell_1.w - cell_2.w);

	RHS_rho.x = RHS_rho.x + RHS.x * d_rho;
	RHS_rho.y = RHS_rho.y + RHS.y  * d_rho;
	RHS_rho.z = RHS_rho.z + RHS.z * d_rho;

	RHS_u.x = RHS_u.x + RHS.x  * d_u;
	RHS_u.y = RHS_u.y + RHS.y  * d_u;
	RHS_u.z = RHS_u.z + RHS.z * d_u;

	RHS_v.x = RHS_v.x + RHS.x  * d_v;
	RHS_v.y = RHS_v.y + RHS.y * d_v;
	RHS_v.z = RHS_v.z + RHS.z * d_v;

	RHS_w.x = RHS_w.x + RHS.x * d_w;
	RHS_w.y = RHS_w.y + RHS.y * d_w;
	RHS_w.z = RHS_w.z + RHS.z* d_w;

}





__global__ void time_integration(int n_cells, int rk, int rk_max_t, double* delta_t_local, double4* soln_t0, double4* soln_t1, double4* soln,
	double* res_rho, double* res_u, double* res_v, double* res_w,double delta_t, int gpu_time_stepping) {



	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	for (int i = index; i < n_cells; i += stride) {
		double f1, f2, f3, f4;
		double alpha[4] = { 1.0 , 0.5, 0.5,1.0 };
	// double beta[4] = { 1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0 };
		 double beta[4] = { 1.0 , 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0 }; // for euler stepping
		if (i < n_cells) {

			double4 t0 = soln_t0[i];
			if (gpu_time_stepping == 1) {
				// update intermediate macroscopic variables for next Runge Kutta Time Step
				f1 = t0.w + res_rho[i] * delta_t_local[i] * alpha[rk];
				f2 = t0.x + (res_u[i]) *delta_t_local[i] * alpha[rk];
				f3 = t0.y + (res_v[i])* delta_t_local[i] * alpha[rk];
				f4 = t0.z + (res_w[i]) * delta_t_local[i] * alpha[rk];
			}
			else {
				f1 = t0.w + res_rho[i] * delta_t * alpha[rk];
				f2 = t0.x + (res_u[i]) * delta_t * alpha[rk];
				f3 = t0.y + (res_v[i])* delta_t * alpha[rk];
				f4 = t0.z + (res_w[i]) * delta_t * alpha[rk];

			}
			

			// change momentum to velocity
			f2 = f2 / f1;
			f3 = f3 / f1;
			f4 = f4 / f1;

			//add contributions to
			double4 t1 = soln_t1[i];
			if (gpu_time_stepping ==1) {
				t1.w = t1.w + delta_t_local[i] * beta[rk] * res_rho[i];
				t1.x = t1.x + delta_t_local[i] * beta[rk] * res_u[i];
				t1.y = t1.y + delta_t_local[i] * beta[rk] * res_v[i];
				t1.z = t1.z + delta_t_local[i] * beta[rk] * res_w[i];
			}
			else {
				t1.w = t1.w + delta_t * beta[rk] * res_rho[i];
				t1.x = t1.x + delta_t * beta[rk] * res_u[i];
				t1.y = t1.y + delta_t * beta[rk] * res_v[i];
				t1.z = t1.z + delta_t * beta[rk] * res_w[i];

			}
			soln_t1[i] = t1; //add update

			double4 solution = soln[i];

			if (rk == (rk_max_t - 1)) {  // assume rk4
				f1 = t1.w;
				f2 = t1.x / t1.w;
				f3 = t1.y / t1.w;
				f4 = t1.z / t1.w;

				solution.w = f1;
				solution.x = f2;
				solution.y = f3;
				solution.z = f4;

			}
			else {
				solution.w = f1;
				solution.x = f2;
				solution.y = f3;
				solution.z = f4;
			}
			soln[i] = solution;


		}

	}



}




__global__ void get_interior_gradients(int n_cells, int* gradient_cells, double4* input_soln, double3* RHS_array,
	double* LHS_xx, double* LHS_xy, double* LHS_xz,
	double* LHS_yx, double* LHS_yy, double* LHS_yz,
	double* LHS_zx, double* LHS_zy, double* LHS_zz,
	double3* grad_rho_arr, double3* grad_u_arr, double3* grad_v_arr, double3* grad_w_arr) {



	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	for (int i = index; i < n_cells; i += stride) {

		if (i < n_cells) {
			int i1;

			double3 grad_u, grad_v, grad_w, grad_rho;

			double3 RHS_u = make_double3(0, 0, 0);
			double3 RHS_v = make_double3(0, 0, 0);
			double3 RHS_w = make_double3(0, 0, 0);
			double3 RHS_rho = make_double3(0, 0, 0);

			//change from CPU code, assume 6 surfaces i.e. hex only.
			// no need for hanging nodes etc. in remainder of project.
			//changed as vectors of vectors is not allowable in CUDA
			for (int i_nb = 0; i_nb < 6; i_nb++) {
				i1 = gradient_cells[i * 6 + i_nb];

				add_LS_contributions(RHS_rho, RHS_u, RHS_v, RHS_w, input_soln, i1, i, i_nb, RHS_array);


			}


			grad_rho.x = LHS_xx[i] * RHS_rho.x + LHS_xy[i] * RHS_rho.y + LHS_xz[i] * RHS_rho.z;
			grad_rho.y = LHS_yx[i] * RHS_rho.x + LHS_yy[i] * RHS_rho.y + LHS_yz[i] * RHS_rho.z;
			grad_rho.z = LHS_zx[i] * RHS_rho.x + LHS_zy[i] * RHS_rho.y + LHS_zz[i] * RHS_rho.z;

			grad_u.x = LHS_xx[i] * RHS_u.x + LHS_xy[i] * RHS_u.y + LHS_xz[i] * RHS_u.z;
			grad_u.y = LHS_yx[i] * RHS_u.x + LHS_yy[i] * RHS_u.y + LHS_yz[i] * RHS_u.z;
			grad_u.z = LHS_zx[i] * RHS_u.x + LHS_zy[i] * RHS_u.y + LHS_zz[i] * RHS_u.z;

			grad_v.x = LHS_xx[i] * RHS_v.x + LHS_xy[i] * RHS_v.y + LHS_xz[i] * RHS_v.z;
			grad_v.y = LHS_yx[i] * RHS_v.x + LHS_yy[i] * RHS_v.y + LHS_yz[i] * RHS_v.z;
			grad_v.z = LHS_zx[i] * RHS_v.x + LHS_zy[i] * RHS_v.y + LHS_zz[i] * RHS_v.z;

			grad_w.x = LHS_xx[i] * RHS_w.x + LHS_xy[i] * RHS_w.y + LHS_xz[i] * RHS_w.z;
			grad_w.y = LHS_yx[i] * RHS_w.x + LHS_yy[i] * RHS_w.y + LHS_yz[i] * RHS_w.z;
			grad_w.z = LHS_zx[i] * RHS_w.x + LHS_zy[i] * RHS_w.y + LHS_zz[i] * RHS_w.z;

			grad_rho_arr[i] = grad_rho;
			grad_u_arr[i] = grad_u;
			grad_v_arr[i] = grad_v;
			grad_w_arr[i] = grad_w;

		}




	}

}

__global__ void add_face_flux_to_cell(int n_faces, double4 * res, double * cell_volume, int* mesh_owner, int* mesh_neighbour, double* res_rho, double* res_u, double* res_v, double* res_w) {

	//loop through cells


	int index = blockIdx.x * blockDim.x + threadIdx.x;;
	int stride = blockDim.x * gridDim.x;



	for (int face = index; face < n_faces; face += stride) {

		if (face < n_faces) {
			int i = mesh_owner[face];
			int neighbour = mesh_neighbour[face];
			double4 cell_flux = res[face];


			//atomic adds due to lack of double in cuda10, use custom function
			//// add density flux to current cell and neighbouring cell
			myatomicAdd(&res_rho[i], -cell_flux.w / cell_volume[i]);
			myatomicAdd(&res_rho[neighbour], cell_flux.w / cell_volume[neighbour]);
			//// add x momentum
			myatomicAdd(&res_u[i], -cell_flux.x / cell_volume[i]);
			myatomicAdd(&res_u[neighbour], cell_flux.x / cell_volume[neighbour]);

			//// add y momentum
			myatomicAdd(&res_v[i], -cell_flux.y / cell_volume[i]);
			myatomicAdd(&res_v[neighbour], cell_flux.y / cell_volume[neighbour]);
			//// add z momentum
			myatomicAdd(&res_w[i], -cell_flux.z / cell_volume[i]);
			myatomicAdd(&res_w[neighbour], cell_flux.z / cell_volume[neighbour]);



		}

	}

}

__global__ void
calc_face_flux_x(int n, double4* input_soln, double* cell_volume, double* surface_area, int* mesh_owner, int* mesh_neighbour, double3* centroid, double3* face_centroid, double3* face_normal,
	double* streaming_dt,
	double3* grad_rho_arr, double3* grad_u_arr, double3* grad_v_arr, double3* grad_w_arr, int n_cells, double pre_conditioned_gamma, double* local_fneq, double visc,
	double* res_rho, double* res_u, double* res_v, double* res_w, double4* res_face,
	int bcs_rho_type[], int bcs_vel_type[], double4 bcs_arr[], double PI, double * visc_factor, int testcase) {

	//loop through cells

	int index = blockIdx.x * blockDim.x + threadIdx.x;;
	int stride = blockDim.x * gridDim.x;



	for (int face = index; face < n; face += stride) {

		if (face < n) {
 
			int i = mesh_owner[face];
			double interface_area = surface_area[face];
			
			int neighbour = mesh_neighbour[face];
			double3 cell_normal = face_normal[face];
			
			double dt = streaming_dt[face]; // dt for the cell interface
			double rho_interface;
			double3 u_interface;
			double4 cell_flux = make_double4(0, 0, 0, 0);

			double cell_1 = centroid[i].x;  //get current cell centre
			double interface_node = face_centroid[face].x;
			double cell_2 = centroid[neighbour].x;

		

			double3 grad_rho_1 = grad_rho_arr[i];
			double3 grad_rho_2 = grad_rho_arr[neighbour];
			double3 grad_u_1 = grad_u_arr[i];
			double3 grad_u_2 = grad_u_arr[neighbour];
			double3 grad_v_1 = grad_v_arr[i];
			double3 grad_v_2 = grad_v_arr[neighbour];
			double3 grad_w_1 = grad_w_arr[i];
			double3 grad_w_2 = grad_w_arr[neighbour];

			double4 owner_soln = input_soln[i];
			double4 neighbour_soln = input_soln[neighbour];

			double u_lattice[15], v_lattice[15], w_lattice[15], rho_lattice[15], feq_lattice[15];

			//populate macro variables
			populate_lattice_macros_uniform_x(u_lattice, v_lattice, w_lattice, rho_lattice, cell_1, cell_2,
				interface_node, i, neighbour, grad_u_1, grad_u_2, grad_v_1, grad_v_2, grad_w_1, grad_w_2, grad_rho_1, grad_rho_2, owner_soln, neighbour_soln, dt);
			populate_feq_lattice_explicit(u_lattice, v_lattice, w_lattice, rho_lattice,
				feq_lattice);



			////get initial feqs
			

			////// get macroscopic values at cell interface
			rho_interface = feq_lattice[0] + feq_lattice[1] + feq_lattice[2] + feq_lattice[3]
				+ feq_lattice[4] + feq_lattice[5] + feq_lattice[6] + feq_lattice[7] + feq_lattice[8]
				+ feq_lattice[9] + feq_lattice[10] + feq_lattice[11] + feq_lattice[12] + feq_lattice[13]
				+ feq_lattice[14];

			u_interface.x = 1 / rho_interface * (feq_lattice[1] - feq_lattice[2]
				+ feq_lattice[7] - feq_lattice[8]
				+ feq_lattice[9] - feq_lattice[10] + feq_lattice[11] - feq_lattice[12] - feq_lattice[13]
				+ feq_lattice[14]);

			u_interface.y = 1 / rho_interface * (feq_lattice[3] - feq_lattice[4]
				+ feq_lattice[7] - feq_lattice[8]
				+ feq_lattice[9] - feq_lattice[10] - feq_lattice[11] + feq_lattice[12] + feq_lattice[13]
				- feq_lattice[14]);

			u_interface.z = 1 / rho_interface * (feq_lattice[5] - feq_lattice[6]
				+ feq_lattice[7] - feq_lattice[8]
				- feq_lattice[9] + feq_lattice[10] + feq_lattice[11] - feq_lattice[12] + feq_lattice[13]
				- feq_lattice[14]);

			//calculate_flux_at_interface(u_interface, dt, pre_conditioned_gamma, rho_interface,
			//	cell_flux, i, feq_lattice, cell_normal, interface_area, local_fneq, visc);
			calculate_flux_at_interface_x(u_interface, dt, pre_conditioned_gamma, rho_interface,
					cell_flux, i, feq_lattice, cell_normal, interface_area, local_fneq, visc,  visc_factor[i], testcase);

			res_face[face] = cell_flux;

		}

	}

}


__global__ void
calc_face_flux_y(int n, double4* input_soln, double* cell_volume, double* surface_area, int* mesh_owner, int* mesh_neighbour, double3* centroid, double3* face_centroid, double3* face_normal,
	double* streaming_dt,
	double3* grad_rho_arr, double3* grad_u_arr, double3* grad_v_arr, double3* grad_w_arr, int n_cells, double pre_conditioned_gamma, double* local_fneq, double visc,
	double* res_rho, double* res_u, double* res_v, double* res_w, double4* res_face,
	int bcs_rho_type[], int bcs_vel_type[], double4 bcs_arr[], double PI,  int start_index, double * visc_factor, int testcase) {

	//loop through cells

	int index = blockIdx.x * blockDim.x + threadIdx.x;;
	int stride = blockDim.x * gridDim.x;



	for (int j = index; j < n; j += stride) {

		if (j < n) {
			int face = start_index + j; // start_index 

			int i = mesh_owner[face];
			double interface_area = surface_area[face];

			int neighbour = mesh_neighbour[face];
			double3 cell_normal = face_normal[face];

			double dt = streaming_dt[face]; // dt for the cell interface
			double rho_interface;
			double3 u_interface;
			double4 cell_flux = make_double4(0, 0, 0, 0);

			double cell_1 = centroid[i].y;  //get current cell centre
			double interface_node = face_centroid[face].y;
			double cell = centroid[neighbour].y;



			double3 grad_rho_1 = grad_rho_arr[i];
			double3 grad_rho = grad_rho_arr[neighbour];
			double3 grad_u_1 = grad_u_arr[i];
			double3 grad_u = grad_u_arr[neighbour];
			double3 grad_v_1 = grad_v_arr[i];
			double3 grad_v = grad_v_arr[neighbour];
			double3 grad_w_1 = grad_w_arr[i];
			double3 grad_w = grad_w_arr[neighbour];

			double4 owner_soln = input_soln[i];
			double4 neighbour_soln = input_soln[neighbour];

			double u_lattice[15], v_lattice[15], w_lattice[15], rho_lattice[15], feq_lattice[15];

			//populate macro variables
			populate_lattice_macros_uniform_y(u_lattice, v_lattice, w_lattice, rho_lattice, cell_1, cell,
				interface_node, i, neighbour, grad_u_1, grad_u, grad_v_1, grad_v, grad_w_1, grad_w, grad_rho_1, grad_rho, owner_soln, neighbour_soln, dt);
			populate_feq_lattice_explicit(u_lattice, v_lattice, w_lattice, rho_lattice,
				feq_lattice);



			////get initial feqs


			////// get macroscopic values at cell interface
			rho_interface = feq_lattice[0] + feq_lattice[1] + feq_lattice[2] + feq_lattice[3]
				+ feq_lattice[4] + feq_lattice[5] + feq_lattice[6] + feq_lattice[7] + feq_lattice[8]
				+ feq_lattice[9] + feq_lattice[10] + feq_lattice[11] + feq_lattice[12] + feq_lattice[13]
				+ feq_lattice[14];

			u_interface.x = 1 / rho_interface * (feq_lattice[1] - feq_lattice[2]
				+ feq_lattice[7] - feq_lattice[8]
				+ feq_lattice[9] - feq_lattice[10] + feq_lattice[11] - feq_lattice[12] - feq_lattice[13]
				+ feq_lattice[14]);

			u_interface.y = 1 / rho_interface * (feq_lattice[3] - feq_lattice[4]
				+ feq_lattice[7] - feq_lattice[8]
				+ feq_lattice[9] - feq_lattice[10] - feq_lattice[11] + feq_lattice[12] + feq_lattice[13]
				- feq_lattice[14]);

			u_interface.z = 1 / rho_interface * (feq_lattice[5] - feq_lattice[6]
				+ feq_lattice[7] - feq_lattice[8]
				- feq_lattice[9] + feq_lattice[10] + feq_lattice[11] - feq_lattice[12] + feq_lattice[13]
				- feq_lattice[14]);

			//calculate_flux_at_interface(u_interface, dt, pre_conditioned_gamma, rho_interface,
			//	cell_flux, i, feq_lattice, cell_normal, interface_area, local_fneq, visc);
			calculate_flux_at_interface_y(u_interface, dt, pre_conditioned_gamma, rho_interface,
				cell_flux, i, feq_lattice, cell_normal, interface_area, local_fneq, visc,visc_factor[i],testcase);

			res_face[face] = cell_flux;

		}

	}

}




__global__ void
calc_face_flux_z(int n, double4* input_soln, double* cell_volume, double* surface_area, int* mesh_owner, int* mesh_neighbour, double3* centroid, double3* face_centroid, double3* face_normal,
	double* streaming_dt,
	double3* grad_rho_arr, double3* grad_u_arr, double3* grad_v_arr, double3* grad_w_arr, int n_cells, double pre_conditioned_gamma, double* local_fneq, double visc,
	double* res_rho, double* res_u, double* res_v, double* res_w, double4* res_face,
	int bcs_rho_type[], int bcs_vel_type[], double4 bcs_arr[], double PI, int start_index, double * visc_factor, int testcase) {

	//loop through cells

	int index = blockIdx.x * blockDim.x + threadIdx.x;;
	int stride = blockDim.x * gridDim.x;



	for (int j = index; j < n; j += stride) {

		if (j < n) {
			int face = start_index + j; // start_index 

			int i = mesh_owner[face];
			double interface_area = surface_area[face];

			int neighbour = mesh_neighbour[face];
			double3 cell_normal = face_normal[face];

			double dt = streaming_dt[face]; // dt for the cell interface
			double rho_interface;
			double3 u_interface;
			double4 cell_flux = make_double4(0, 0, 0, 0);

			double cell_1 = centroid[i].z;  //get current cell centre
			double interface_node = face_centroid[face].z;
			double cell = centroid[neighbour].z;



			double3 grad_rho_1 = grad_rho_arr[i];
			double3 grad_rho = grad_rho_arr[neighbour];
			double3 grad_u_1 = grad_u_arr[i];
			double3 grad_u = grad_u_arr[neighbour];
			double3 grad_v_1 = grad_v_arr[i];
			double3 grad_v = grad_v_arr[neighbour];
			double3 grad_w_1 = grad_w_arr[i];
			double3 grad_w = grad_w_arr[neighbour];

			double4 owner_soln = input_soln[i];
			double4 neighbour_soln = input_soln[neighbour];

			double u_lattice[15], v_lattice[15], w_lattice[15], rho_lattice[15], feq_lattice[15];

			//populate macro variables
			populate_lattice_macros_uniform_z(u_lattice, v_lattice, w_lattice, rho_lattice, cell_1, cell,
				interface_node, i, neighbour, grad_u_1, grad_u, grad_v_1, grad_v, grad_w_1, grad_w, grad_rho_1, grad_rho, owner_soln, neighbour_soln, dt);
			populate_feq_lattice_explicit(u_lattice, v_lattice, w_lattice, rho_lattice,
				feq_lattice);



			////get initial feqs


			////// get macroscopic values at cell interface
			rho_interface = feq_lattice[0] + feq_lattice[1] + feq_lattice[2] + feq_lattice[3]
				+ feq_lattice[4] + feq_lattice[5] + feq_lattice[6] + feq_lattice[7] + feq_lattice[8]
				+ feq_lattice[9] + feq_lattice[10] + feq_lattice[11] + feq_lattice[12] + feq_lattice[13]
				+ feq_lattice[14];

			u_interface.x = 1 / rho_interface * (feq_lattice[1] - feq_lattice[2]
				+ feq_lattice[7] - feq_lattice[8]
				+ feq_lattice[9] - feq_lattice[10] + feq_lattice[11] - feq_lattice[12] - feq_lattice[13]
				+ feq_lattice[14]);

			u_interface.y = 1 / rho_interface * (feq_lattice[3] - feq_lattice[4]
				+ feq_lattice[7] - feq_lattice[8]
				+ feq_lattice[9] - feq_lattice[10] - feq_lattice[11] + feq_lattice[12] + feq_lattice[13]
				- feq_lattice[14]);

			u_interface.z = 1 / rho_interface * (feq_lattice[5] - feq_lattice[6]
				+ feq_lattice[7] - feq_lattice[8]
				- feq_lattice[9] + feq_lattice[10] + feq_lattice[11] - feq_lattice[12] + feq_lattice[13]
				- feq_lattice[14]);

			//calculate_flux_at_interface(u_interface, dt, pre_conditioned_gamma, rho_interface,
			//	cell_flux, i, feq_lattice, cell_normal, interface_area, local_fneq, visc);
			calculate_flux_at_interface_z(u_interface, dt, pre_conditioned_gamma, rho_interface,
				cell_flux, i, feq_lattice, cell_normal, interface_area, local_fneq, visc,visc_factor[i],testcase);

			res_face[face] = cell_flux;

		}

	}

}




__global__ void
calc_face_flux(int n, double4* input_soln, double* cell_volume, double* surface_area, int* mesh_owner, int* mesh_neighbour, double3* centroid, double3* face_centroid, double3* face_normal,
	double* streaming_dt,
	double3* grad_rho_arr, double3* grad_u_arr, double3* grad_v_arr, double3* grad_w_arr, int n_cells, double pre_conditioned_gamma, double* local_fneq, double visc,
	double* res_rho, double* res_u, double* res_v, double* res_w, double4* res_face,
	int bcs_rho_type[], int bcs_vel_type[], double4 bcs_arr[], double PI, int testcase, double * visc_factor, bool preconditioning) {

	//loop through cells

	int index = blockIdx.x * blockDim.x + threadIdx.x;;
	int stride = blockDim.x * gridDim.x;

	__shared__ double d_1;
	d_1 = sin(PI / 14); // distance comparisons to see if lattice nodes are within 5 degrees of the shared face
	__shared__ double d_173;
	d_173 = sin(PI / 14)*sqrt(3.0); // distance comparisons to see if lattice nodes are within 5 degrees of the shared face



	for (int face = index; face < n; face += stride) {

		if (face < n) {
			

			int i = mesh_owner[face];
			double interface_area = surface_area[face];
			double3 cell_1 = centroid[i]; //get current cell centre
			double3 interface_node = face_centroid[face];
			int neighbour = mesh_neighbour[face];
			double3 cell_normal = face_normal[face];
			double3 cell_2 = centroid[neighbour];
			double4 interface_macrovariables = make_double4(0, 0, 0, 0);
			double dt = streaming_dt[face]; // dt for the cell interface
			
			double3 grad_rho_1 = grad_rho_arr[i];
			double3 grad_rho_2 = grad_rho_arr[neighbour];
			double3 grad_u_1 = grad_u_arr[i];
			double3 grad_u_2 = grad_u_arr[neighbour];
			double3 grad_v_1 = grad_v_arr[i];
			double3 grad_v_2 = grad_v_arr[neighbour];
			double3 grad_w_1 = grad_w_arr[i];
			double3 grad_w_2 = grad_w_arr[neighbour];

			double4 owner_soln = input_soln[i];
			double4 neighbour_soln = input_soln[neighbour];

			double rho_interface;
			double3 u_interface;

			double4 cell_flux = make_double4(0, 0, 0, 0);

			double u_lattice[15], v_lattice[15], w_lattice[15], rho_lattice[15], feq_lattice[15];



			//populate macro variables
			populate_lattice_macros(u_lattice, v_lattice, w_lattice, rho_lattice, cell_1, cell_2,
				interface_node, i, neighbour, grad_u_1, grad_u_2, grad_v_1, grad_v_2, grad_w_1, grad_w_2, grad_rho_1, grad_rho_2, owner_soln, neighbour_soln, n_cells, cell_normal, dt,
				bcs_rho_type, bcs_vel_type, bcs_arr, d_1, d_173);
			//populate macro variables


			//get initial feqs
			if (preconditioning) {
				for (int k = 0; k < 15; k++) {
					populate_feq(u_lattice, v_lattice, w_lattice, rho_lattice,
						feq_lattice, k, pre_conditioned_gamma);
				}

			}
			else {
				populate_feq_lattice_explicit(u_lattice, v_lattice, w_lattice, rho_lattice,
					feq_lattice);

			}
			
			
				
			
			// get macroscopic values at cell interface
			rho_interface = feq_lattice[0] + feq_lattice[1] + feq_lattice[2] + feq_lattice[3]
				+ feq_lattice[4] + feq_lattice[5] + feq_lattice[6] + feq_lattice[7] + feq_lattice[8]
				+ feq_lattice[9] + feq_lattice[10] + feq_lattice[11] + feq_lattice[12] + feq_lattice[13]
				+ feq_lattice[14];

			u_interface.x = 1 / rho_interface * (feq_lattice[1] - feq_lattice[2]
				+ feq_lattice[7] - feq_lattice[8]
				+ feq_lattice[9] - feq_lattice[10] + feq_lattice[11] - feq_lattice[12] - feq_lattice[13]
				+ feq_lattice[14]);

			u_interface.y = 1 / rho_interface * (feq_lattice[3] - feq_lattice[4]
				+ feq_lattice[7] - feq_lattice[8]
				+ feq_lattice[9] - feq_lattice[10] - feq_lattice[11] + feq_lattice[12] + feq_lattice[13]
				- feq_lattice[14]);

			u_interface.z = 1 / rho_interface * (feq_lattice[5] - feq_lattice[6]
				+ feq_lattice[7] - feq_lattice[8]
				- feq_lattice[9] + feq_lattice[10] + feq_lattice[11] - feq_lattice[12] + feq_lattice[13]
				- feq_lattice[14]);
			

				calculate_flux_at_interface(u_interface, dt, pre_conditioned_gamma, rho_interface,
					cell_flux, i, feq_lattice, cell_normal, interface_area, local_fneq, visc, testcase, visc_factor[i]);

			

			res_face[face] = cell_flux;
					

		}

	}

}



__global__ void
calc_face_flux_bc(int n, double4* input_soln, double* cell_volume, double* surface_area, int* mesh_owner, int* mesh_neighbour, double3* centroid, double3* face_centroid, double3* face_normal,
	double* streaming_dt,
	double3* grad_rho_arr, double3* grad_u_arr, double3* grad_v_arr, double3* grad_w_arr, int n_cells, double pre_conditioned_gamma, double* local_fneq, double visc,
	double* res_rho, double* res_u, double* res_v, double* res_w, double4* res_face,
	int bcs_rho_type[], int bcs_vel_type[], double4 bcs_arr[], double PI, int start_index, int testcase, bool preconditioning ) {

	//loop through cells

	int index = blockIdx.x * blockDim.x + threadIdx.x;;
	int stride = blockDim.x * gridDim.x;

	__shared__ double d_1;
	d_1 = sin(PI / 14); // distance comparisons to see if lattice nodes are within 5 degrees of the shared face
	__shared__ double d_173;
	d_173 = sin(PI / 14)*sqrt(3.0); // distance comparisons to see if lattice nodes are within 5 degrees of the shared face



	for (int j = index; j < n; j += stride) {

		if (j < n) {
			int face = j + start_index;

			int i = mesh_owner[face];
			double interface_area = surface_area[face];
			double3 cell_1 = centroid[i]; //get current cell centre
			double3 interface_node = face_centroid[face];
			int neighbour = mesh_neighbour[face];
			double3 cell_normal = face_normal[face];
			double3 cell = centroid[neighbour];
			double4 interface_macrovariables = make_double4(0, 0, 0, 0);
			double dt = streaming_dt[face]; // dt for the cell interface

			double3 grad_rho_1 = grad_rho_arr[i];
			double3 grad_rho = grad_rho_arr[neighbour];
			double3 grad_u_1 = grad_u_arr[i];
			double3 grad_u = grad_u_arr[neighbour];
			double3 grad_v_1 = grad_v_arr[i];
			double3 grad_v = grad_v_arr[neighbour];
			double3 grad_w_1 = grad_w_arr[i];
			double3 grad_w = grad_w_arr[neighbour];

			double4 owner_soln = input_soln[i];
			double4 neighbour_soln = input_soln[neighbour];

			double rho_interface;
			double3 u_interface;

			double4 cell_flux = make_double4(0, 0, 0, 0);

			double u_lattice[15], v_lattice[15], w_lattice[15], rho_lattice[15], feq_lattice[15];



			//populate macro variables
/*		populate_lattice_macros_uniform(u_lattice, v_lattice, w_lattice, rho_lattice, cell_1, cell,
			interface_node, i, neighbour, grad_u_1, grad_u, grad_v_1, grad_v , grad_w_1, grad_w, grad_rho_1, grad_rho, owner_soln, neighbour_soln, n_cells,cell_normal,dt,
			 bcs_rho_type,  bcs_vel_type, bcs_arr);*/

			populate_lattice_macros(u_lattice, v_lattice, w_lattice, rho_lattice, cell_1, cell,
				interface_node, i, neighbour, grad_u_1, grad_u, grad_v_1, grad_v, grad_w_1, grad_w, grad_rho_1, grad_rho, owner_soln, neighbour_soln, n_cells, cell_normal, dt,
				bcs_rho_type, bcs_vel_type, bcs_arr, d_1, d_173);
			//populate macro variables


			//get initial feqs
			if (preconditioning) {
				for (int k = 0; k < 15; k++) {
					populate_feq(u_lattice, v_lattice, w_lattice, rho_lattice,
						feq_lattice, k, pre_conditioned_gamma);
				}

			}
			else {
				populate_feq_lattice_explicit(u_lattice, v_lattice, w_lattice, rho_lattice,
					feq_lattice);

			}


			// get macroscopic values at cell interface
			rho_interface = feq_lattice[0] + feq_lattice[1] + feq_lattice[2] + feq_lattice[3]
				+ feq_lattice[4] + feq_lattice[5] + feq_lattice[6] + feq_lattice[7] + feq_lattice[8]
				+ feq_lattice[9] + feq_lattice[10] + feq_lattice[11] + feq_lattice[12] + feq_lattice[13]
				+ feq_lattice[14];

			u_interface.x = 1 / rho_interface * (feq_lattice[1] - feq_lattice[2]
				+ feq_lattice[7] - feq_lattice[8]
				+ feq_lattice[9] - feq_lattice[10] + feq_lattice[11] - feq_lattice[12] - feq_lattice[13]
				+ feq_lattice[14]);

			u_interface.y = 1 / rho_interface * (feq_lattice[3] - feq_lattice[4]
				+ feq_lattice[7] - feq_lattice[8]
				+ feq_lattice[9] - feq_lattice[10] - feq_lattice[11] + feq_lattice[12] + feq_lattice[13]
				- feq_lattice[14]);

			u_interface.z = 1 / rho_interface * (feq_lattice[5] - feq_lattice[6]
				+ feq_lattice[7] - feq_lattice[8]
				- feq_lattice[9] + feq_lattice[10] + feq_lattice[11] - feq_lattice[12] + feq_lattice[13]
				- feq_lattice[14]);

			

			calculate_flux_at_interface(u_interface, dt, pre_conditioned_gamma, rho_interface,
				cell_flux, i, feq_lattice, cell_normal, interface_area, local_fneq, visc, testcase, 1.0);

			


			res_face[face] = cell_flux;

			


		}

	}

}


__global__ void calc_total_residual(int n_cells, double4 res_tot, double* res_rho, double* res_u, double* res_v, double* res_w) {

	//loop through cells

	int index = blockIdx.x * blockDim.x + threadIdx.x;;
	int stride = blockDim.x * gridDim.x;

	for (int i = index; i < n_cells; i += stride) {
		if (i < n_cells) {
			myatomicAdd(&res_tot.w, res_rho[i] * res_rho[i]);
			myatomicAdd(&res_tot.x, res_u[i] * res_u[i]);
			myatomicAdd(&res_tot.y, res_v[i] * res_v[i]);
			myatomicAdd(&res_tot.z, res_w[i] * res_w[i]);

		}
	}
}



__global__ void check_error(int n_cells, double4* soln) {

	//loop through cells

	int index = blockIdx.x * blockDim.x + threadIdx.x;;
	int stride = blockDim.x * gridDim.x;

	for (int i = index; i < n_cells; i += stride) {
		if (i < n_cells) {
			double4 tmp = soln[i];

			if (isnan(tmp.w) || isnan(tmp.x)) {
				printf("nan failure");

				/*printf << "nan failure" << endl;
				cout << t << endl;
				cout << i << endl;*/

				asm("trap;");
			}

			if (tmp.w > 1000) {
				printf("rho failure");
				/*	cout << "rho failure" << endl;
					cout << t << endl;
					cout << i << endl;*/

				asm("trap;");
			}

		}
	}
}





__global__ void get_cfl_device(int n, double4* input_soln, double* cell_volume, double* delta_t_local, double3* cfl_areas, double factor,
	double max_velocity, double pre_conditioned_gamma, double visc, int gpu_time_stepping) {

	//loop through cells

	int index = blockIdx.x * blockDim.x + threadIdx.x;;
	int stride = blockDim.x * gridDim.x;

	for (int i = index; i < n; i += stride) {
		double effective_speed_of_sound;
		double area_x_eigen, visc_eigen;
		double visc_constant;
		visc_constant = 4;
		double  temp;
		double3 area;
		double4 soln;
		area = cfl_areas[i];
		soln = input_soln[i];		//effective_speed_of_sound = 1/sqrt(3);
		effective_speed_of_sound = max_velocity * sqrt(1 - pre_conditioned_gamma +
			pow(pre_conditioned_gamma / sqrt(3.0) / max_velocity, 2));

		if (i < n) {


			// eigen values as per Zhaoli guo(2004) - preconditioning

			  //estimation of spectral radii s per Jiri Blasek: CFD Principles and Application Determination of Max time Step
			/*area_x_eigen = cell_volume[i];
			area_x_eigen = 0;*/
			/*area_x_eigen = (soln.x + effective_speed_of_sound)*area.x
				+ (soln.y + effective_speed_of_sound)*area.y
				+ (soln.z + effective_speed_of_sound)*area.z;*/

			area_x_eigen = (max_velocity + effective_speed_of_sound)*area.x
				+ (max_velocity + effective_speed_of_sound)*area.y
				+ (max_velocity + effective_speed_of_sound)*area.z;

			area_x_eigen = area_x_eigen / pre_conditioned_gamma;

			//reducing preconditioning increases viscous flux - increases eigenvalue
			visc_eigen = 2 * visc / pre_conditioned_gamma / soln.w / cell_volume[i];
			visc_eigen = visc_eigen * (area.x * area.x + area.y * area.y + area.z * area.z);

			area_x_eigen = area_x_eigen + visc_constant * visc_eigen;

			// use smallest time step allowed
			temp = factor * cell_volume[i] / area_x_eigen;
			

			//fully local time stepping
			delta_t_local[i] = temp;

			//if (gpu_time_stepping == 1 || gpu_time_stepping == 4) {
			//	delta_t_local[i] = temp;

			//}
			//else { //constant user defined time step
			//	delta_t_local[i] = factor;

			//}

		}
	}

	return;
}

//owner is always west and normal is west to east
__device__ void populate_lattice_macros_uniform_x(double u_lattice[], double v_lattice[],
	double w_lattice[], double rho_lattice[], double cell_1, double cell_2,
	double interface_node, int i, int neighbour,
	double3 grad_u_1, double3 grad_u_2, double3 grad_v_1, double3  grad_v_2, double3  grad_w_1, double3  grad_w_2, double3  grad_rho_1, double3  grad_rho_2,
	double4 owner_soln, double4 neighbour_soln, double dt) {

	//double3 temp1, temp2;
	double temp1, temp2;
	double rho_i, rho_nb, u_i, u_nb, v_i, v_nb, w_i, w_nb;
	///   case 0: // center node

	temp1 = cell_1 - interface_node;
	temp2 = cell_2 - interface_node;

	double3 grad_rho, grad_u, grad_v, grad_w;
	
	rho_i = owner_soln.w - temp1 * grad_rho_1.x;
	rho_nb = neighbour_soln.w - (temp2* grad_rho_2.x);
	rho_lattice[0] = (rho_i + rho_nb)*0.5;

	u_i = owner_soln.x - temp1 * grad_u_1.x;
	u_nb = neighbour_soln.x - (temp2* grad_u_2.x);
	u_lattice[0] = (u_i + u_nb)*0.5;

	v_i = owner_soln.y  - temp1 * grad_v_1.x;
	v_nb = neighbour_soln.y - (temp2* grad_v_2.x);
	v_lattice[0] = (v_i + v_nb)*0.5;

	w_i = owner_soln.z - temp1 * grad_w_1.x;
	w_nb = neighbour_soln.z - (temp2* grad_w_2.x);
	w_lattice[0] = (w_i + w_nb)*0.5;
		
	double a, b, c;
	c = cell_2 - cell_1;
	a = fabs(temp2 / c);
	b = fabs(temp1 / c);
	
	grad_u.x = (grad_u_1.x *a + grad_u_2.x *b);
	grad_u.y = (grad_u_1.y *a + grad_u_2.y *b);
	grad_u.z = (grad_u_1.z *a + grad_u_2.z *b);

	grad_v.x = (grad_v_1.x *a + grad_v_2.x *b);
	grad_v.y = (grad_v_1.y *a + grad_v_2.y *b);
	grad_v.z = (grad_v_1.z *a + grad_v_2.z *b);

	grad_w.x = (grad_w_1.x *a + grad_w_2.x *b);
	grad_w.y = (grad_w_1.y *a + grad_w_2.y *b);
	grad_w.z = (grad_w_1.z *a + grad_w_2.z *b);

	grad_rho.x = (grad_rho_1.x *a + grad_rho_2.x *b);
	grad_rho.y = (grad_rho_1.y *a + grad_rho_2.y *b);
	grad_rho.z = (grad_rho_1.z *a + grad_rho_2.z *b);


	   
	///  case 1:west_node


	rho_lattice[1] = rho_i - grad_rho_1.x* dt;
	u_lattice[1] = u_i - grad_u_1.x* dt;
	v_lattice[1] = v_i - grad_v_1.x *dt;
	w_lattice[1] = w_i - grad_w_1.x *dt;



	///  case 2: // east_node

	rho_lattice[2] = rho_nb + grad_rho_2.x* dt;
	u_lattice[2] = u_nb + grad_u_2.x* dt;
	v_lattice[2] = v_nb + grad_v_2.x *dt;
	w_lattice[2] = w_nb + grad_w_2.x *dt;
	   
	///   case 3: // bottom node

	rho_lattice[3] = rho_lattice[0] - grad_rho.y* dt;
	u_lattice[3] = u_lattice[0] - grad_u.y* dt;
	v_lattice[3] = v_lattice[0] - grad_v.y *dt;
	w_lattice[3] = w_lattice[0] - grad_w.y *dt;
	
	///   case 4: // top node
	
	rho_lattice[4] = rho_lattice[0] + grad_rho.y* dt;
	u_lattice[4] = u_lattice[0] + grad_u.y* dt;
	v_lattice[4] = v_lattice[0] + grad_v.y *dt;
	w_lattice[4] = w_lattice[0] + grad_w.y *dt;

	///   case 5: // back node

	rho_lattice[5] = rho_lattice[0] - grad_rho.z* dt;
	u_lattice[5] = u_lattice[0] - grad_u.z* dt;
	v_lattice[5] = v_lattice[0] - grad_v.z *dt;
	w_lattice[5] = w_lattice[0] - grad_w.z *dt;

	///   case 6: // front node

	rho_lattice[6] = rho_lattice[0] + grad_rho.z* dt;
	u_lattice[6] = u_lattice[0] + grad_u.z* dt;
	v_lattice[6] = v_lattice[0] + grad_v.z *dt;
	w_lattice[6] = w_lattice[0] + grad_w.z *dt;



	/// case 7: back bottom west

	rho_lattice[7] = rho_lattice[1] 
		- grad_rho_1.y* dt
		- grad_rho_1.z* dt;
	u_lattice[7] = u_lattice[1] 
		- grad_u_1.y* dt
		- grad_u_1.z* dt;
	v_lattice[7] = v_lattice[1] 
		- grad_v_1.y* dt
		- grad_v_1.z* dt;
	w_lattice[7] = w_lattice[1] 
		- grad_w_1.y* dt
		- grad_w_1.z* dt;

	/// case 9: front bottom west

	rho_lattice[9] = rho_lattice[1]
		- grad_rho_1.y* dt
		+ grad_rho_1.z* dt;
	u_lattice[9] = u_lattice[1]
		- grad_u_1.y* dt
		+ grad_u_1.z* dt;
	v_lattice[9] = v_lattice[1]
		- grad_v_1.y* dt
		+ grad_v_1.z* dt;
	w_lattice[9] = w_lattice[1]
		- grad_w_1.y* dt
		+ grad_w_1.z* dt;


	///  case 11: back top west

	rho_lattice[11] = rho_lattice[1]
		+ grad_rho_1.y* dt
		- grad_rho_1.z* dt;
	u_lattice[11] = u_lattice[1]
		+ grad_u_1.y* dt
		- grad_u_1.z* dt;
	v_lattice[11] = v_lattice[1]
		+ grad_v_1.y* dt
		- grad_v_1.z* dt;
	w_lattice[11] = w_lattice[1]
		+ grad_w_1.y* dt
		- grad_w_1.z* dt;


	/// case 14: front top west

	rho_lattice[14] = rho_lattice[1]
		+ grad_rho_1.y* dt
		+ grad_rho_1.z* dt;
	u_lattice[14] = u_lattice[1]
		+ grad_u_1.y* dt
		+ grad_u_1.z* dt;
	v_lattice[14] = v_lattice[1]
		+ grad_v_1.y* dt
		+ grad_v_1.z* dt;
	w_lattice[14] = w_lattice[1]
		+ grad_w_1.y* dt
		+ grad_w_1.z* dt;
	   


	/// case 8: front top east

	rho_lattice[8] = rho_lattice[2] 
		+ grad_rho_2.y* dt
		+ grad_rho_2.z* dt;
	u_lattice[8] = u_lattice[2] 
		+ grad_u_2.y* dt
		+ grad_u_2.z* dt;
	v_lattice[8] = v_lattice[2] 
		+ grad_v_2.y* dt
		+ grad_v_2.z* dt;
	w_lattice[8] = w_lattice[2] 
		+ grad_w_2.y* dt
		+ grad_w_2.z* dt;



	/// case 10 Back Top East

	rho_lattice[10] = rho_lattice[2]
		+ grad_rho_2.y* dt
		- grad_rho_2.z* dt;
	u_lattice[10] = u_lattice[2]
		+ grad_u_2.y* dt
		- grad_u_2.z* dt;
	v_lattice[10] = v_lattice[2]
		+ grad_v_2.y* dt
		- grad_v_2.z* dt;
	w_lattice[10] = w_lattice[2]
		+ grad_w_2.y* dt
		- grad_w_2.z* dt;


	/// case 12 Front Bottom East

	rho_lattice[12] = rho_lattice[2]
		- grad_rho_2.y* dt
		+ grad_rho_2.z* dt;
	u_lattice[12] = u_lattice[2]
		- grad_u_2.y* dt
		+ grad_u_2.z* dt;
	v_lattice[12] = v_lattice[2]
		- grad_v_2.y* dt
		+ grad_v_2.z* dt;
	w_lattice[12] = w_lattice[2]
		- grad_w_2.y* dt
		+ grad_w_2.z* dt;
	   


	/// case 13 Back Bottom East
	rho_lattice[13] = rho_lattice[2]
		- grad_rho_2.y* dt
		- grad_rho_2.z* dt;
	u_lattice[13] = u_lattice[2]
		- grad_u_2.y* dt
		- grad_u_2.z* dt;
	v_lattice[13] = v_lattice[2]
		- grad_v_2.y* dt
		- grad_v_2.z* dt;
	w_lattice[13] = w_lattice[2]
		- grad_w_2.y* dt
		- grad_w_2.z* dt;



}



//owner is always west and normal is west to east
__device__ void populate_lattice_macros_uniform_y(double u_lattice[], double v_lattice[],
	double w_lattice[], double rho_lattice[], double cell_1, double cell_2,
	double interface_node, int i, int neighbour,
	double3 grad_u_1, double3 grad_u_2, double3 grad_v_1, double3  grad_v_2, double3  grad_w_1, double3  grad_w_2, double3  grad_rho_1, double3  grad_rho_2,
	double4 owner_soln, double4 neighbour_soln, double dt) {

	//double3 temp1, temp2;
	double temp1, temp2;
	double rho_i, rho_nb, u_i, u_nb, v_i, v_nb, w_i, w_nb;
	///   case 0: // center node

	temp1 = cell_1 - interface_node;
	temp2 = cell_2 - interface_node;

	double3 grad_rho, grad_u, grad_v, grad_w;
	
	rho_i = owner_soln.w - temp1 * grad_rho_1.y;
	rho_nb = neighbour_soln.w - (temp2* grad_rho_2.y);
	rho_lattice[0] = (rho_i + rho_nb)*0.5;

	u_i = owner_soln.x - temp1 * grad_u_1.y;
	u_nb = neighbour_soln.x - (temp2* grad_u_2.y);
	u_lattice[0] = (u_i + u_nb)*0.5;

	v_i = owner_soln.y  - temp1 * grad_v_1.y;
	v_nb = neighbour_soln.y - (temp2* grad_v_2.y);
	v_lattice[0] = (v_i + v_nb)*0.5;

	w_i = owner_soln.z - temp1 * grad_w_1.y;
	w_nb = neighbour_soln.z - (temp2* grad_w_2.y);
	w_lattice[0] = (w_i + w_nb)*0.5;
		
	double a, b, c;
	c = cell_2 - cell_1;
	a = fabs(temp2 / c);
	b = fabs(temp1 / c);
	
	grad_u.x = (grad_u_1.x *a + grad_u_2.x *b);
	grad_u.y = (grad_u_1.y *a + grad_u_2.y *b);
	grad_u.z = (grad_u_1.z *a + grad_u_2.z *b);

	grad_v.x = (grad_v_1.x *a + grad_v_2.x *b);
	grad_v.y = (grad_v_1.y *a + grad_v_2.y *b);
	grad_v.z = (grad_v_1.z *a + grad_v_2.z *b);

	grad_w.x = (grad_w_1.x *a + grad_w_2.x *b);
	grad_w.y = (grad_w_1.y *a + grad_w_2.y *b);
	grad_w.z = (grad_w_1.z *a + grad_w_2.z *b);

	grad_rho.x = (grad_rho_1.x *a + grad_rho_2.x *b);
	grad_rho.y = (grad_rho_1.y *a + grad_rho_2.y *b);
	grad_rho.z = (grad_rho_1.z *a + grad_rho_2.z *b);


	   
	///  case 1:west_node


	rho_lattice[1] = rho_lattice[0] - grad_rho.x* dt;
	u_lattice[1] = u_lattice[0] - grad_u.x* dt;
	v_lattice[1] = v_lattice[0] - grad_v.x *dt;
	w_lattice[1] = w_lattice[0] - grad_w.x *dt;



	///  case 2: // east_node

	rho_lattice[2] = rho_lattice[0] + grad_rho.x* dt;
	u_lattice[2] = u_lattice[0] + grad_u.x* dt;
	v_lattice[2] = v_lattice[0] + grad_v.x *dt;
	w_lattice[2] = w_lattice[0] + grad_w.x *dt;
	   
	///   case 3: // bottom node

	rho_lattice[3] = rho_i - grad_rho_1.y* dt;
	u_lattice[3] = u_i - grad_u_1.y* dt;
	v_lattice[3] = v_i - grad_v_1.y *dt;
	w_lattice[3] = w_i - grad_w_1.y *dt;
	
	///   case 4: // top node
	
	rho_lattice[4] = rho_nb + grad_rho_2.y* dt;
	u_lattice[4] = u_nb + grad_u_2.y* dt;
	v_lattice[4] = v_nb + grad_v_2.y *dt;
	w_lattice[4] = w_nb + grad_w_2.y *dt;

	///   case 5: // back node

	rho_lattice[5] = rho_lattice[0] - grad_rho.z* dt;
	u_lattice[5] = u_lattice[0] - grad_u.z* dt;
	v_lattice[5] = v_lattice[0] - grad_v.z *dt;
	w_lattice[5] = w_lattice[0] - grad_w.z *dt;

	///   case 6: // front node

	rho_lattice[6] = rho_lattice[0] + grad_rho.z* dt;
	u_lattice[6] = u_lattice[0] + grad_u.z* dt;
	v_lattice[6] = v_lattice[0] + grad_v.z *dt;
	w_lattice[6] = w_lattice[0] + grad_w.z *dt;



	/// case 7: back bottom west

	rho_lattice[7] = rho_lattice[3] 
		- grad_rho_1.x* dt
		- grad_rho_1.z* dt;
	u_lattice[7] = u_lattice[3] 
		- grad_u_1.x* dt
		- grad_u_1.z* dt;
	v_lattice[7] = v_lattice[3] 
		- grad_v_1.x* dt
		- grad_v_1.z* dt;
	w_lattice[7] = w_lattice[3] 
		- grad_w_1.x* dt
		- grad_w_1.z* dt;

	/// case 9: front bottom west

	rho_lattice[9] = rho_lattice[3]
		- grad_rho_1.x* dt
		+ grad_rho_1.z* dt;
	u_lattice[9] = u_lattice[3]
		- grad_u_1.x* dt
		+ grad_u_1.z* dt;
	v_lattice[9] = v_lattice[3]
		- grad_v_1.x* dt
		+ grad_v_1.z* dt;
	w_lattice[9] = w_lattice[3]
		- grad_w_1.x* dt
		+ grad_w_1.z* dt;



	/// case 12 Front Bottom East

	rho_lattice[12] = rho_lattice[3]
		+ grad_rho_1.x* dt
		+ grad_rho_1.z* dt;
	u_lattice[12] = u_lattice[3]
		+ grad_u_1.x* dt
		+ grad_u_1.z* dt;
	v_lattice[12] = v_lattice[3]
		+ grad_v_1.x* dt
		+ grad_v_1.z* dt;
	w_lattice[12] = w_lattice[3]
		+ grad_w_1.x* dt
		+ grad_w_1.z* dt;
	   


	/// case 13 Back Bottom East

	rho_lattice[13] = rho_lattice[3]
		+ grad_rho_1.x* dt
		- grad_rho_1.z* dt;
	u_lattice[13] = u_lattice[3]
		+ grad_u_1.x* dt
		- grad_u_1.z* dt;
	v_lattice[13] = v_lattice[3]
		+ grad_v_1.x* dt
		- grad_v_1.z* dt;
	w_lattice[13] = w_lattice[3]
		+ grad_w_1.x* dt
		- grad_w_1.z* dt;

	

	///  case 11: back top west
	
	rho_lattice[11] = rho_lattice[4]
		- grad_rho_2.x* dt
		- grad_rho_2.z* dt;
	u_lattice[11] = u_lattice[4]
		- grad_u_2.x* dt
		- grad_u_2.z* dt;
	v_lattice[11] = v_lattice[4]
		- grad_v_2.x* dt
		- grad_v_2.z* dt;
	w_lattice[11] = w_lattice[4]
		- grad_w_2.x* dt
		- grad_w_2.z* dt;


	/// case 14: front top west

	rho_lattice[14] = rho_lattice[4]
		- grad_rho_2.x* dt
		+ grad_rho_2.z* dt;
	u_lattice[14] = u_lattice[4]
		- grad_u_2.x* dt
		+ grad_u_2.z* dt;
	v_lattice[14] = v_lattice[4]
		- grad_v_2.x* dt
		+ grad_v_2.z* dt;
	w_lattice[14] = w_lattice[4]
		- grad_w_2.x* dt
		+ grad_w_2.z* dt;

	   
	/// case 8: front top east


	rho_lattice[8] = rho_lattice[4]
		+ grad_rho_2.x* dt
		+ grad_rho_2.z* dt;
	u_lattice[8] = u_lattice[4]
		+ grad_u_2.x* dt
		+ grad_u_2.z* dt;
	v_lattice[8] = v_lattice[4]
		+ grad_v_2.x* dt
		+ grad_v_2.z* dt;
	w_lattice[8] = w_lattice[4]
		+ grad_w_2.x* dt
		+ grad_w_2.z* dt;

	
	/// case 10 Back Top East

	rho_lattice[10] = rho_lattice[4]
		+ grad_rho_2.x* dt
		- grad_rho_2.z* dt;
	u_lattice[10] = u_lattice[4]
		+ grad_u_2.x* dt
		- grad_u_2.z* dt;
	v_lattice[10] = v_lattice[4]
		+ grad_v_2.x* dt
		- grad_v_2.z* dt;
	w_lattice[10] = w_lattice[4]
		+ grad_w_2.x* dt
		- grad_w_2.z* dt;


}


//owner is always west and normal is west to east
__device__ void populate_lattice_macros_uniform_z(double u_lattice[], double v_lattice[],
	double w_lattice[], double rho_lattice[], double cell_1, double cell_2,
	double interface_node, int i, int neighbour,
	double3 grad_u_1, double3 grad_u_2, double3 grad_v_1, double3  grad_v_2, double3  grad_w_1, double3  grad_w_2, double3  grad_rho_1, double3  grad_rho_2,
	double4 owner_soln, double4 neighbour_soln, double dt) {

	//double3 temp1, temp2;
	double temp1, temp2;
	double rho_i, rho_nb, u_i, u_nb, v_i, v_nb, w_i, w_nb;
	///   case 0: // center node

	temp1 = cell_1 - interface_node;
	temp2 = cell_2 - interface_node;

	double3 grad_rho, grad_u, grad_v, grad_w;

	rho_i = owner_soln.w - temp1 * grad_rho_1.z;
	rho_nb = neighbour_soln.w - (temp2* grad_rho_2.z);
	rho_lattice[0] = (rho_i + rho_nb)*0.5;

	u_i = owner_soln.x - temp1 * grad_u_1.z;
	u_nb = neighbour_soln.x - (temp2* grad_u_2.z);
	u_lattice[0] = (u_i + u_nb)*0.5;

	v_i = owner_soln.y - temp1 * grad_v_1.z;
	v_nb = neighbour_soln.y - (temp2* grad_v_2.z);
	v_lattice[0] = (v_i + v_nb)*0.5;

	w_i = owner_soln.z - temp1 * grad_w_1.z;
	w_nb = neighbour_soln.z - (temp2* grad_w_2.z);
	w_lattice[0] = (w_i + w_nb)*0.5;

	double a, b, c;
	c = cell_2 - cell_1;
	a = fabs(temp2 / c);
	b = fabs(temp1 / c);

	grad_u.x = (grad_u_1.x *a + grad_u_2.x *b);
	grad_u.y = (grad_u_1.y *a + grad_u_2.y *b);
	grad_u.z = (grad_u_1.z *a + grad_u_2.z *b);

	grad_v.x = (grad_v_1.x *a + grad_v_2.x *b);
	grad_v.y = (grad_v_1.y *a + grad_v_2.y *b);
	grad_v.z = (grad_v_1.z *a + grad_v_2.z *b);

	grad_w.x = (grad_w_1.x *a + grad_w_2.x *b);
	grad_w.y = (grad_w_1.y *a + grad_w_2.y *b);
	grad_w.z = (grad_w_1.z *a + grad_w_2.z *b);

	grad_rho.x = (grad_rho_1.x *a + grad_rho_2.x *b);
	grad_rho.y = (grad_rho_1.y *a + grad_rho_2.y *b);
	grad_rho.z = (grad_rho_1.z *a + grad_rho_2.z *b);



	///  case 1:west_node


	rho_lattice[1] = rho_lattice[0] - grad_rho.x* dt;
	u_lattice[1] = u_lattice[0] - grad_u.x* dt;
	v_lattice[1] = v_lattice[0] - grad_v.x *dt;
	w_lattice[1] = w_lattice[0] - grad_w.x *dt;
	
	///  case 2: // east_node

	rho_lattice[2] = rho_lattice[0] + grad_rho.x* dt;
	u_lattice[2] = u_lattice[0] + grad_u.x* dt;
	v_lattice[2] = v_lattice[0] + grad_v.x *dt;
	w_lattice[2] = w_lattice[0] + grad_w.x *dt;

	///   case 3: // bottom node

	rho_lattice[3] = rho_lattice[0] - grad_rho.y* dt;
	u_lattice[3] = u_lattice[0] - grad_u.y* dt;
	v_lattice[3] = v_lattice[0] - grad_v.y *dt;
	w_lattice[3] = w_lattice[0] - grad_w.y *dt;

	///   case 4: // top node

	rho_lattice[4] = rho_lattice[0] + grad_rho.y* dt;
	u_lattice[4] = u_lattice[0] + grad_u.y* dt;
	v_lattice[4] = v_lattice[0] + grad_v.y *dt;
	w_lattice[4] = w_lattice[0] + grad_w.y *dt;

	///   case 5: // back node

	rho_lattice[5] = rho_i - grad_rho_1.z* dt;
	u_lattice[5] = u_i - grad_u_1.z* dt;
	v_lattice[5] = v_i - grad_v_1.z *dt;
	w_lattice[5] = w_i - grad_w_1.z *dt;

	///   case 6: // front node

	rho_lattice[6] = rho_nb + grad_rho_2.z* dt;
	u_lattice[6] = u_nb + grad_u_2.z* dt;
	v_lattice[6] = v_nb + grad_v_2.z *dt;
	w_lattice[6] = w_nb + grad_w_2.z *dt;


	/// case 7: back bottom west

	rho_lattice[7] = rho_lattice[5]
		- grad_rho_1.x* dt
		- grad_rho_1.y* dt;
	u_lattice[7] = u_lattice[5]
		- grad_u_1.x* dt
		- grad_u_1.y* dt;
	v_lattice[7] = v_lattice[5]
		- grad_v_1.x* dt
		- grad_v_1.y* dt;
	w_lattice[7] = w_lattice[5]
		- grad_w_1.x* dt
		- grad_w_1.y* dt;


	/// case 13 Back Bottom East

	rho_lattice[13] = rho_lattice[5]
		+ grad_rho_1.x* dt
		- grad_rho_1.y* dt;
	u_lattice[13] = u_lattice[5]
		+ grad_u_1.x* dt
		- grad_u_1.y* dt;
	v_lattice[13] = v_lattice[5]
		+ grad_v_1.x* dt
		- grad_v_1.y* dt;
	w_lattice[13] = w_lattice[5]
		+ grad_w_1.x* dt
		- grad_w_1.y* dt;


	/// case 10 Back Top East

	rho_lattice[10] = rho_lattice[5]
		+ grad_rho_1.x* dt
		+ grad_rho_1.y* dt;
	u_lattice[10] = u_lattice[5]
		+ grad_u_1.x* dt
		+ grad_u_1.y* dt;
	v_lattice[10] = v_lattice[5]
		+ grad_v_1.x* dt
		+ grad_v_1.y* dt;
	w_lattice[10] = w_lattice[5]
		+ grad_w_1.x* dt
		+ grad_w_1.y* dt;




	///  case 11: back top west

	rho_lattice[11] = rho_lattice[5]
		- grad_rho_1.x* dt
		+ grad_rho_1.y* dt;
	u_lattice[11] = u_lattice[5]
		- grad_u_1.x* dt
		+ grad_u_1.y* dt;
	v_lattice[11] = v_lattice[5]
		- grad_v_1.x* dt
		+ grad_v_1.y* dt;
	w_lattice[11] = w_lattice[5]
		- grad_w_1.x* dt
		+ grad_w_1.y* dt;


	/// case 9: front bottom west

	rho_lattice[9] = rho_lattice[6]
		- grad_rho_2.x* dt
		- grad_rho_2.y* dt;
	u_lattice[9] = u_lattice[6]
		- grad_u_2.x* dt
		- grad_u_2.y* dt;
	v_lattice[9] = v_lattice[6]
		- grad_v_2.x* dt
		- grad_v_2.y* dt;
	w_lattice[9] = w_lattice[6]
		- grad_w_2.x* dt
		- grad_w_2.y* dt;

	
	/// case 12 Front Bottom East


	rho_lattice[12] = rho_lattice[6]
		+ grad_rho_2.x* dt
		- grad_rho_2.y* dt;
	u_lattice[12] = u_lattice[6]
		+ grad_u_2.x* dt
		- grad_u_2.y* dt;
	v_lattice[12] = v_lattice[6]
		+ grad_v_2.x* dt
		- grad_v_2.y* dt;
	w_lattice[12] = w_lattice[6]
		+ grad_w_2.x* dt
		- grad_w_2.y* dt;
	   	 

	/// case 14: front top west

	rho_lattice[14] = rho_lattice[6]
		- grad_rho_2.x* dt
		+ grad_rho_2.y* dt;
	u_lattice[14] = u_lattice[6]
		- grad_u_2.x* dt
		+ grad_u_2.y* dt;
	v_lattice[14] = v_lattice[6]
		- grad_v_2.x* dt
		+ grad_v_2.y* dt;
	w_lattice[14] = w_lattice[6]
		- grad_w_2.x* dt
		+ grad_w_2.y* dt;



	/// case 8: front top east

	rho_lattice[8] = rho_lattice[6]
		+ grad_rho_2.x* dt
		+ grad_rho_2.y* dt;
	u_lattice[8] = u_lattice[6]
		+ grad_u_2.x* dt
		+ grad_u_2.y* dt;
	v_lattice[8] = v_lattice[6]
		+ grad_v_2.x* dt
		+ grad_v_2.y* dt;
	w_lattice[8] = w_lattice[6]
		+ grad_w_2.x* dt
		+ grad_w_2.y* dt;


}



__device__ void calculate_flux_at_interface_x(double3 u_interface, double dt, double pre_conditioned_gamma,
	double rho_interface, double4 &cell_flux, int i,
	double* feq_lattice, double3 cell_normal, double interface_area, double* local_fneq, double visc, double visc_factor, int testcase) {

	double uu2, vv2, ww2, u2v2w2, uu, vv, ww, uv, uw, vw, fneq_tau;
	double feq_interface[15];
	double4 x_flux;

	uu2 = u_interface.x * u_interface.x ;
	vv2 = u_interface.y * u_interface.y ;
	ww2 = u_interface.z *u_interface.z ;

	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	uu = u_interface.x;
	vv = u_interface.y;
	ww = u_interface.z;

	uv = uu * vv*9.0 ;
	uw = uu * ww*9.0 ;
	vw = vv * ww *9.0 ;


	if (testcase == 8) {
		fneq_tau = 0.0;
	}
	else {
		fneq_tau = (visc* visc_factor * 3 / dt * pre_conditioned_gamma);
	}
	// local_fneq[i] = fneq_tau;


	feq_interface[1] = lattice_weight[1] * rho_interface*
		(1.0 + 3.0*uu + 4.5*uu2 - u2v2w2);
	feq_interface[1] = feq_interface[1]
		- fneq_tau * (feq_interface[1] - feq_lattice[1]);

	feq_interface[2] = lattice_weight[2] * rho_interface*
		(1.0 - 3.0*uu + 4.5*uu2 - u2v2w2);
	feq_interface[2] = feq_interface[2]
		- fneq_tau * (feq_interface[2] - feq_lattice[2]);

	/*feq_interface[3] = lattice_weight[3] * rho_interface*
		(1.0 + 3.0*vv + 4.5*vv2 - u2v2w2);
	feq_interface[3] = feq_interface[3]
		- fneq_tau * (feq_interface[3] - feq_lattice[3]);

	feq_interface[4] = lattice_weight[4] * rho_interface*
		(1.0 - 3.0*vv + 4.5*vv2 - u2v2w2);
	feq_interface[4] = feq_interface[4]
		- fneq_tau * (feq_interface[4] - feq_lattice[4]);

	feq_interface[5] = lattice_weight[5] * rho_interface*
		(1.0 + 3.0*ww + 4.5*ww2 - u2v2w2);
	feq_interface[5] = feq_interface[5]
		- fneq_tau * (feq_interface[5] - feq_lattice[5]);

	feq_interface[6] = lattice_weight[6] * rho_interface*
		(1.0 - 3.0*ww + 4.5*ww2 - u2v2w2);
	feq_interface[6] = feq_interface[6]
		- fneq_tau * (feq_interface[6] - feq_lattice[6]);*/

	feq_interface[7] = lattice_weight[7] * rho_interface*
		(1.0 + 3.0*uu + 3.0*vv + 3.0*ww + uv + uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[7] = feq_interface[7]
		- fneq_tau * (feq_interface[7] - feq_lattice[7]);

	feq_interface[8] = lattice_weight[8] * rho_interface*
		(1.0 - 3.0*uu - 3.0*vv - 3.0*ww + uv + uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[8] = feq_interface[8]
		- fneq_tau * (feq_interface[8] - feq_lattice[8]);

	feq_interface[9] = lattice_weight[9] * rho_interface*
		(1.0 + 3.0*uu + 3.0*vv - 3.0*ww + uv - uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[9] = feq_interface[9]
		- fneq_tau * (feq_interface[9] - feq_lattice[9]);

	feq_interface[10] = lattice_weight[10] * rho_interface*
		(1.0 - 3.0*uu - 3.0*vv + 3.0*ww + uv - uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[10] = feq_interface[10]
		- fneq_tau * (feq_interface[10] - feq_lattice[10]);

	feq_interface[11] = lattice_weight[11] * rho_interface*
		(1.0 + 3.0*uu - 3.0*vv + 3.0*ww - uv + uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[11] = feq_interface[11]
		- fneq_tau * (feq_interface[11] - feq_lattice[11]);

	feq_interface[12] = lattice_weight[12] * rho_interface*
		(1.0 - 3.0*uu + 3.0*vv - 3.0*ww - uv + uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[12] = feq_interface[12]
		- fneq_tau * (feq_interface[12] - feq_lattice[12]);

	feq_interface[13] = lattice_weight[13] * rho_interface*
		(1.0 - 3.0*uu + 3.0*vv + 3.0*ww - uv - uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[13] = feq_interface[13]
		- fneq_tau * (feq_interface[13] - feq_lattice[13]);

	feq_interface[14] = lattice_weight[14] * rho_interface*
		(1.0 + 3.0*uu - 3.0*vv - 3.0*ww - uv - uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[14] = feq_interface[14]
		- fneq_tau * (feq_interface[14] - feq_lattice[14]);


	x_flux.w = (feq_interface[1] - feq_interface[2]
		+ feq_interface[7] - feq_interface[8]
		+ feq_interface[9] - feq_interface[10] + feq_interface[11]
		- feq_interface[12] - feq_interface[13]
		+ feq_interface[14]);

	x_flux.x =
		(feq_interface[1] + feq_interface[2]
			+ feq_interface[7] + feq_interface[8]
			+ feq_interface[9] + feq_interface[10] + feq_interface[11]
			+ feq_interface[12] + feq_interface[13]
			+ feq_interface[14]);

	x_flux.y =
		feq_interface[7] + feq_interface[8]
		+ feq_interface[9] + feq_interface[10] - feq_interface[11] - feq_interface[12]
		- feq_interface[13]
		- feq_interface[14];

	x_flux.z =
		feq_interface[7] + feq_interface[8]
		- feq_interface[9] - feq_interface[10] + feq_interface[11] + feq_interface[12]
		- feq_interface[13]
		- feq_interface[14];
	   	 
	cell_flux.w = (x_flux.w )*interface_area;
	cell_flux.x = (x_flux.x 
		)*interface_area;

	cell_flux.y = (x_flux.y 
		)*interface_area;


	cell_flux.z = (x_flux.z 
		)*interface_area;
}

__device__ void calculate_flux_at_interface_y(double3 u_interface, double dt, double pre_conditioned_gamma,
	double rho_interface, double4 &cell_flux, int i,
	double* feq_lattice, double3 cell_normal, double interface_area, double* local_fneq, double visc, double visc_factor,int testcase) {

	double uu2, vv2, ww2, u2v2w2, uu, vv, ww, uv, uw, vw, fneq_tau;
	double feq_interface[15];
	double4  y_flux;

	uu2 = u_interface.x * u_interface.x;
	vv2 = u_interface.y * u_interface.y;
	ww2 = u_interface.z *u_interface.z;

	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	uu = u_interface.x;
	vv = u_interface.y;
	ww = u_interface.z;

	uv = uu * vv*9.0;
	uw = uu * ww*9.0;
	vw = vv * ww *9.0;

	if (testcase == 8) {
		fneq_tau = 0.0;
	}
	else {
		fneq_tau = (visc* visc_factor * 3 / dt * pre_conditioned_gamma);
	}

	
	// local_fneq[i] = fneq_tau;


	//feq_interface[1] = lattice_weight[1] * rho_interface*
	//	(1.0 + 3.0*uu + 4.5*uu2 - u2v2w2);
	//feq_interface[1] = feq_interface[1]
	//	- fneq_tau * (feq_interface[1] - feq_lattice[1]);

	//feq_interface[2] = lattice_weight[2] * rho_interface*
	//	(1.0 - 3.0*uu + 4.5*uu2 - u2v2w2);
	//feq_interface[2] = feq_interface[2]
	//	- fneq_tau * (feq_interface[2] - feq_lattice[2]);

	feq_interface[3] = lattice_weight[3] * rho_interface*
		(1.0 + 3.0*vv + 4.5*vv2 - u2v2w2);
	feq_interface[3] = feq_interface[3]
		- fneq_tau * (feq_interface[3] - feq_lattice[3]);

	feq_interface[4] = lattice_weight[4] * rho_interface*
		(1.0 - 3.0*vv + 4.5*vv2 - u2v2w2);
	feq_interface[4] = feq_interface[4]
		- fneq_tau * (feq_interface[4] - feq_lattice[4]);

	//feq_interface[5] = lattice_weight[5] * rho_interface*
	//	(1.0 + 3.0*ww + 4.5*ww2 - u2v2w2);
	//feq_interface[5] = feq_interface[5]
	//	- fneq_tau * (feq_interface[5] - feq_lattice[5]);

	//feq_interface[6] = lattice_weight[6] * rho_interface*
	//	(1.0 - 3.0*ww + 4.5*ww2 - u2v2w2);
	//feq_interface[6] = feq_interface[6]
	//	- fneq_tau * (feq_interface[6] - feq_lattice[6]);

	feq_interface[7] = lattice_weight[7] * rho_interface*
		(1.0 + 3.0*uu + 3.0*vv + 3.0*ww + uv + uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[7] = feq_interface[7]
		- fneq_tau * (feq_interface[7] - feq_lattice[7]);

	feq_interface[8] = lattice_weight[8] * rho_interface*
		(1.0 - 3.0*uu - 3.0*vv - 3.0*ww + uv + uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[8] = feq_interface[8]
		- fneq_tau * (feq_interface[8] - feq_lattice[8]);

	feq_interface[9] = lattice_weight[9] * rho_interface*
		(1.0 + 3.0*uu + 3.0*vv - 3.0*ww + uv - uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[9] = feq_interface[9]
		- fneq_tau * (feq_interface[9] - feq_lattice[9]);

	feq_interface[10] = lattice_weight[10] * rho_interface*
		(1.0 - 3.0*uu - 3.0*vv + 3.0*ww + uv - uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[10] = feq_interface[10]
		- fneq_tau * (feq_interface[10] - feq_lattice[10]);

	feq_interface[11] = lattice_weight[11] * rho_interface*
		(1.0 + 3.0*uu - 3.0*vv + 3.0*ww - uv + uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[11] = feq_interface[11]
		- fneq_tau * (feq_interface[11] - feq_lattice[11]);

	feq_interface[12] = lattice_weight[12] * rho_interface*
		(1.0 - 3.0*uu + 3.0*vv - 3.0*ww - uv + uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[12] = feq_interface[12]
		- fneq_tau * (feq_interface[12] - feq_lattice[12]);

	feq_interface[13] = lattice_weight[13] * rho_interface*
		(1.0 - 3.0*uu + 3.0*vv + 3.0*ww - uv - uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[13] = feq_interface[13]
		- fneq_tau * (feq_interface[13] - feq_lattice[13]);

	feq_interface[14] = lattice_weight[14] * rho_interface*
		(1.0 + 3.0*uu - 3.0*vv - 3.0*ww - uv - uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[14] = feq_interface[14]
		- fneq_tau * (feq_interface[14] - feq_lattice[14]);


	y_flux.w = (feq_interface[3] - feq_interface[4]
		+ feq_interface[7] - feq_interface[8]
		+ feq_interface[9] - feq_interface[10] - feq_interface[11]
		+ feq_interface[12] + feq_interface[13]
		- feq_interface[14]);


	y_flux.x = feq_interface[7] + feq_interface[8]
		+ feq_interface[9] + feq_interface[10] - feq_interface[11] - feq_interface[12]
		- feq_interface[13]
		- feq_interface[14];

	y_flux.y = (feq_interface[3] + feq_interface[4]
		+ feq_interface[7] + feq_interface[8]
		+ feq_interface[9] + feq_interface[10] + feq_interface[11]
		+ feq_interface[12] + feq_interface[13]
		+ feq_interface[14]);

	y_flux.z = feq_interface[7] + feq_interface[8]
		- feq_interface[9] - feq_interface[10] - feq_interface[11] - feq_interface[12]
		+ feq_interface[13]
		+ feq_interface[14];

	cell_flux.w = (y_flux.w)*interface_area;
	cell_flux.x = (y_flux.x)*interface_area;

	cell_flux.y = (y_flux.y)*interface_area;


	cell_flux.z = (y_flux.z)*interface_area;





}


__device__ void calculate_flux_at_interface_z(double3 u_interface, double dt, double pre_conditioned_gamma,
	double rho_interface, double4 &cell_flux, int i,
	double* feq_lattice, double3 cell_normal, double interface_area, double* local_fneq, double visc, double visc_factor, int testcase) {

	double uu2, vv2, ww2, u2v2w2, uu, vv, ww, uv, uw, vw, fneq_tau;
	double feq_interface[15];
	double4  z_flux;

	uu2 = u_interface.x * u_interface.x;
	vv2 = u_interface.y * u_interface.y;
	ww2 = u_interface.z *u_interface.z;

	u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

	uu = u_interface.x;
	vv = u_interface.y;
	ww = u_interface.z;

	uv = uu * vv*9.0;
	uw = uu * ww*9.0;
	vw = vv * ww *9.0;


	if (testcase == 8) {
		fneq_tau = 0.0;
	}
	else {
		fneq_tau = (visc* visc_factor * 3 / dt * pre_conditioned_gamma);
	}

	// local_fneq[i] = fneq_tau;


	//feq_interface[1] = lattice_weight[1] * rho_interface*
	//	(1.0 + 3.0*uu + 4.5*uu2 - u2v2w2);
	//feq_interface[1] = feq_interface[1]
	//	- fneq_tau * (feq_interface[1] - feq_lattice[1]);

	//feq_interface[2] = lattice_weight[2] * rho_interface*
	//	(1.0 - 3.0*uu + 4.5*uu2 - u2v2w2);
	//feq_interface[2] = feq_interface[2]
	//	- fneq_tau * (feq_interface[2] - feq_lattice[2]);

	//feq_interface[3] = lattice_weight[3] * rho_interface*
	//	(1.0 + 3.0*vv + 4.5*vv2 - u2v2w2);
	//feq_interface[3] = feq_interface[3]
	//	- fneq_tau * (feq_interface[3] - feq_lattice[3]);

	//feq_interface[4] = lattice_weight[4] * rho_interface*
	//	(1.0 - 3.0*vv + 4.5*vv2 - u2v2w2);
	//feq_interface[4] = feq_interface[4]
	//	- fneq_tau * (feq_interface[4] - feq_lattice[4]);

	feq_interface[5] = lattice_weight[5] * rho_interface*
		(1.0 + 3.0*ww + 4.5*ww2 - u2v2w2);
	feq_interface[5] = feq_interface[5]
		- fneq_tau * (feq_interface[5] - feq_lattice[5]);

	feq_interface[6] = lattice_weight[6] * rho_interface*
		(1.0 - 3.0*ww + 4.5*ww2 - u2v2w2);
	feq_interface[6] = feq_interface[6]
		- fneq_tau * (feq_interface[6] - feq_lattice[6]);

	feq_interface[7] = lattice_weight[7] * rho_interface*
		(1.0 + 3.0*uu + 3.0*vv + 3.0*ww + uv + uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[7] = feq_interface[7]
		- fneq_tau * (feq_interface[7] - feq_lattice[7]);

	feq_interface[8] = lattice_weight[8] * rho_interface*
		(1.0 - 3.0*uu - 3.0*vv - 3.0*ww + uv + uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[8] = feq_interface[8]
		- fneq_tau * (feq_interface[8] - feq_lattice[8]);

	feq_interface[9] = lattice_weight[9] * rho_interface*
		(1.0 + 3.0*uu + 3.0*vv - 3.0*ww + uv - uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[9] = feq_interface[9]
		- fneq_tau * (feq_interface[9] - feq_lattice[9]);

	feq_interface[10] = lattice_weight[10] * rho_interface*
		(1.0 - 3.0*uu - 3.0*vv + 3.0*ww + uv - uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[10] = feq_interface[10]
		- fneq_tau * (feq_interface[10] - feq_lattice[10]);

	feq_interface[11] = lattice_weight[11] * rho_interface*
		(1.0 + 3.0*uu - 3.0*vv + 3.0*ww - uv + uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[11] = feq_interface[11]
		- fneq_tau * (feq_interface[11] - feq_lattice[11]);

	feq_interface[12] = lattice_weight[12] * rho_interface*
		(1.0 - 3.0*uu + 3.0*vv - 3.0*ww - uv + uw - vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[12] = feq_interface[12]
		- fneq_tau * (feq_interface[12] - feq_lattice[12]);

	feq_interface[13] = lattice_weight[13] * rho_interface*
		(1.0 - 3.0*uu + 3.0*vv + 3.0*ww - uv - uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[13] = feq_interface[13]
		- fneq_tau * (feq_interface[13] - feq_lattice[13]);

	feq_interface[14] = lattice_weight[14] * rho_interface*
		(1.0 + 3.0*uu - 3.0*vv - 3.0*ww - uv - uw + vw +
			4.5*uu2 + 4.5*vv2 + 4.5*ww2 - u2v2w2);
	feq_interface[14] = feq_interface[14]
		- fneq_tau * (feq_interface[14] - feq_lattice[14]);


	z_flux.w = (feq_interface[5] - feq_interface[6]
		+ feq_interface[7] - feq_interface[8]
		- feq_interface[9] + feq_interface[10] + feq_interface[11]
		- feq_interface[12] + feq_interface[13]
		- feq_interface[14]);
	z_flux.x = feq_interface[7] + feq_interface[8]
		- feq_interface[9] - feq_interface[10] + feq_interface[11] + feq_interface[12]
		- feq_interface[13]
		- feq_interface[14];
	z_flux.y = feq_interface[7] + feq_interface[8]
		- feq_interface[9] - feq_interface[10] - feq_interface[11] - feq_interface[12]
		+ feq_interface[13]
		+ feq_interface[14];
	z_flux.z = (feq_interface[5] + feq_interface[6]
		+ feq_interface[7] + feq_interface[8]
		+ feq_interface[9] + feq_interface[10] + feq_interface[11]
		+ feq_interface[12] + feq_interface[13]
		+ feq_interface[14]);





	cell_flux.w = (z_flux.w)*interface_area;
	cell_flux.x = (z_flux.x)*interface_area;

	cell_flux.y = (z_flux.y)*interface_area;


	cell_flux.z = (z_flux.z)*interface_area;





}