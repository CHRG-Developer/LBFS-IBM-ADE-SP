#include <stdlib.h>
#include <math.h>
#include "post_processing.h"
#include <iostream>
#include "gradients.h"
#include <cstdio>
#include <fstream>
#include <sstream>

post_processing::post_processing()
{
    //ctor
}

post_processing::~post_processing()
{
    //dtor

     delete [] vorticity;
    vorticity = NULL;
    delete [] streamfunction;
    streamfunction = NULL;

}
post_processing::post_processing(int _total_nodes)
{
    //ctor
    total_nodes = _total_nodes;
     //rho = (double*) malloc (sizeof(double)*(total_nodes));
     vorticity = new double [total_nodes +1];
        if (vorticity==NULL) exit (1);
     streamfunction = new double [total_nodes+1];
        if (streamfunction==NULL) exit (1);
        Initialise();



}

void post_processing::Initialise() {

    std::fill_n(vorticity, total_nodes , 0.00);
    std::fill_n(streamfunction, total_nodes, 0.0);

    drag_coef = 0.0;
    lift_coef = 0.0;

}

void post_processing::cylinder_post_processing( unstructured_mesh &mesh, global_variables &globals,
                gradients &grads, Boundary_Conditions &bcs, Solution &soln,domain_geometry &geom , Solution &wall_shear_stress){



        //calculate lift and drag forces
        std::string output_location;
    std::string filename;
    std::ofstream drag_text ;
    std::string drag_file;
    output_location = globals.output_file;
    filename = globals.simulation_name;
    drag_file = output_location + "/drag.txt";

/// Generic Load Case Input

     drag_text.open(drag_file.c_str(), std::ios::out);


	 drag_text << "pressure" << " , " << "visc" << " , " << "n_i" << " , " << "n_j" << " , " << "n_k" <<
		 " , " << "u.x" << " , " << "u.y" << " , " << "area" << " , " << "wall_shear_stress.get_rho(wall)"
		 << " , " << "wall_shear_stress.get_u(wall)" << " , " << "wall_shear_stress.get_v(wall)" << " , " << "shear_stress"
		 << " , " << "mesh.get_centroid_x(owner)" << " , " << "mesh.get_centroid_y(owner)" << std::endl;


        double lift, drag;

        lift = 0.0;
        drag = 0.0;

        int face, owner,bc;

        double pressure, visc;

        double n_i,n_j,n_k,area,s_i,s_j,s_k,shear_stress;
		double traction_x, traction_y;

        vector_var u,v, vel,  normal,vel_parrallel, centroid_bc, centroid_nb;

        double normal_grad;
		int wall;
		wall = 0;

    for (int i = 0; i < mesh.get_num_bc(); i++){

        if( mesh.bc_types[i] == "wall"){

            face = i + mesh.get_n_neighbours();
            owner = mesh.get_mesh_owner(face);
            bc = mesh.get_mesh_neighbour(face);

            pressure = soln.get_rho(owner) /3.0 * globals.pre_conditioned_gamma;

            visc = globals.visc * globals.pre_conditioned_gamma;

            u.set_equal(grads.get_u(bc));
            v.set_equal(grads.get_v(bc));

            n_i = mesh.get_face_i(face);
            n_j = mesh.get_face_j(face);
            n_k = mesh.get_face_k(face);
            area = mesh.get_face_area(face);

            drag = drag + ((  pressure - visc* u.x)*n_i
                    - visc * u.y *n_j)* area;

            lift = lift + ( visc*v.x*n_i - (pressure - visc*v.y)*n_j )*area;
			// traction force is momentum without pressure contributions
			// Tx = mom_x - p
			//Ty = mom_y +p
			traction_x = wall_shear_stress.get_u(wall) - wall_shear_stress.get_rho(wall) / 3 * globals.pre_conditioned_gamma*area;
			traction_y = wall_shear_stress.get_v(wall) + wall_shear_stress.get_rho(wall) / 3 * globals.pre_conditioned_gamma*area;
			//parrallel vectors - two choices and both give the same zero point
			s_i = n_j;
			s_j = -n_i;

			//shear stress in wall is dot(s, traction_force)
			shear_stress = s_i * traction_x + s_j * traction_y;


            drag_text << pressure << " , " << visc << " , " << n_i << " , " << n_j << " , " << n_k <<
                                " , " << u.x  << " , " << u.y << " , " << area << " , " << wall_shear_stress.get_rho(wall) 
				<< " , " << wall_shear_stress.get_u(wall) << " , " << wall_shear_stress.get_v(wall) << " , " << shear_stress 
				<< " , " << mesh.get_centroid_x(owner) << " , " << mesh.get_centroid_y(owner) << std::endl;


            //calculate parrallel velocity p. 605 moukallad
            vel.set_equal ( soln.get_u(owner), soln.get_v(owner), soln.get_w(owner));
            normal.set_equal(n_i,n_j,n_k);


            vel_parrallel.set_equal(vel.x - vel.Dot_Product(normal) * n_i,
                                    vel.y - vel.Dot_Product(normal) *n_j,
                                    vel.z - vel.Dot_Product(normal) *n_k);

            centroid_bc.set_equal(mesh.get_centroid_x(bc), mesh.get_centroid_y(bc), mesh.get_centroid_z(bc));
            centroid_nb.set_equal(mesh.get_centroid_x(owner),mesh.get_centroid_y(owner),mesh.get_centroid_z(owner));

            centroid_nb.subtract(centroid_bc);
            normal_grad = vel_parrallel.Magnitude()/ centroid_nb.Magnitude();

            separation_grad.push_back(normal_grad);
            separation_angle.push_back(atan( n_j/n_i) *360/2/globals.PI);
            separation_i.push_back(n_i);
            separation_j.push_back(n_j);
            separation_x.push_back(centroid_bc.x);
            separation_y.push_back(centroid_bc.y);

			wall = wall + 1;


        }

    }

    drag = drag;
    lift = lift;

    drag_coef = 2*drag/globals.ref_rho / pow(globals.max_velocity,2)/geom.Y /geom.Z;
    lift_coef = 2*lift/globals.ref_rho / pow(globals.max_velocity,2)/geom.X /geom.Z;

    drag_text.close();

    }



	void post_processing::cylinder_post_processing(unstructured_mesh &mesh, global_variables &globals,
		gradients &grads, Boundary_Conditions &bcs, Solution &soln, domain_geometry &geom, double4 * res_face) {



		//calculate lift and drag forces
		std::string output_location;
		std::string filename;
		std::ofstream drag_text;
		std::string drag_file;
		output_location = globals.output_file;
		filename = globals.simulation_name;
		drag_file = output_location + "/drag.txt";

		/// Generic Load Case Input

		drag_text.open(drag_file.c_str(), std::ios::out);


		drag_text << "pressure" << " , " << "visc" << " , " << "n_i" << " , " << "n_j" << " , " << "n_k" <<
			" , " << "u.x" << " , " << "u.y" << " , " << "area" << " , " << "wall_shear_stress.get_rho(wall)"
			<< " , " << "wall_shear_stress.get_u(wall)" << " , " << "wall_shear_stress.get_v(wall)" << " , " << "shear_stress"
			<< " , " << "mesh.get_centroid_x(owner)" << " , " << "mesh.get_centroid_y(owner)" << std::endl;


		double lift, drag;

		lift = 0.0;
		drag = 0.0;

		int face, owner, bc;

		double pressure, visc;

		double n_i, n_j, n_k, area, s_i, s_j, s_k, shear_stress;
		double traction_x, traction_y;
		double4 wall_res;

		vector_var u, v, vel, normal, vel_parrallel, centroid_bc, centroid_nb;

		double normal_grad;
		int wall;
		wall = 0;

		for (int i = 0; i < mesh.get_num_bc(); i++) {

			if (mesh.bc_types[i] == "wall") {

				face = i + mesh.get_n_neighbours();
				owner = mesh.get_mesh_owner(face);
				bc = mesh.get_mesh_neighbour(face);

				pressure = soln.get_rho(owner) / 3.0 * globals.pre_conditioned_gamma;

				visc = globals.visc * globals.pre_conditioned_gamma;

				u.set_equal(grads.get_u(bc));
				v.set_equal(grads.get_v(bc));

				n_i = mesh.get_face_i(face);
				n_j = mesh.get_face_j(face);
				n_k = mesh.get_face_k(face);
				area = mesh.get_face_area(face);

				drag = drag + ((pressure - visc * u.x)*n_i
					- visc * u.y *n_j)* area;

				lift = lift + (visc*v.x*n_i - (pressure - visc * v.y)*n_j)*area;

				wall_res = res_face[face];

				// traction force is momentum without pressure contributions
				// Tx = mom_x - p
				//Ty = mom_y +p
			//	traction_x = wall_shear_stress.get_u(wall) - wall_shear_stress.get_rho(wall) / 3 * globals.pre_conditioned_gamma*area;
			//	traction_y = wall_shear_stress.get_v(wall) + wall_shear_stress.get_rho(wall) / 3 * globals.pre_conditioned_gamma*area;
				//parrallel vectors - two choices and both give the same zero point
				s_i = n_j;
				s_j = -n_i;

				//shear stress in wall is dot(s, traction_force)
				shear_stress = s_i * traction_x + s_j * traction_y;


				drag_text << pressure << " , " << visc << " , " << n_i << " , " << n_j << " , " << n_k <<
					" , " << u.x << " , " << u.y << " , " << area << " , " << wall_res.w
					<< " , " << wall_res.x << " , " << wall_res.y << " , " << shear_stress
					<< " , " << mesh.get_centroid_x(owner) << " , " << mesh.get_centroid_y(owner) << std::endl;


				//calculate parrallel velocity p. 605 moukallad
				vel.set_equal(soln.get_u(owner), soln.get_v(owner), soln.get_w(owner));
				normal.set_equal(n_i, n_j, n_k);


				vel_parrallel.set_equal(vel.x - vel.Dot_Product(normal) * n_i,
					vel.y - vel.Dot_Product(normal) *n_j,
					vel.z - vel.Dot_Product(normal) *n_k);

				centroid_bc.set_equal(mesh.get_centroid_x(bc), mesh.get_centroid_y(bc), mesh.get_centroid_z(bc));
				centroid_nb.set_equal(mesh.get_centroid_x(owner), mesh.get_centroid_y(owner), mesh.get_centroid_z(owner));

				centroid_nb.subtract(centroid_bc);
				normal_grad = vel_parrallel.Magnitude() / centroid_nb.Magnitude();

				separation_grad.push_back(normal_grad);
				separation_angle.push_back(atan(n_j / n_i) * 360 / 2 / globals.PI);
				separation_i.push_back(n_i);
				separation_j.push_back(n_j);
				separation_x.push_back(centroid_bc.x);
				separation_y.push_back(centroid_bc.y);

				wall = wall + 1;


			}

		}

		drag = drag;
		lift = lift;

		drag_coef = 2 * drag / globals.ref_rho / pow(globals.max_velocity, 2) / geom.Y / geom.Z;
		lift_coef = 2 * lift / globals.ref_rho / pow(globals.max_velocity, 2) / geom.X / geom.Z;

		drag_text.close();

	}

void post_processing::calc_vorticity(Solution &x_grad, Solution & y_grad){

    for(int i =0; i < total_nodes; i++){
        vorticity[i] = x_grad.get_v(i) - y_grad.get_u(i);

    }

}

void post_processing::calc_streamfunction(Mesh &mesh, global_variables &globals,
            Boundary_Conditions &bcs){


    double t,w;

    t= cos(globals.PI/mesh.get_num_x()) + cos(globals.PI/mesh.get_num_y());

    w = (8 - sqrt(64-16*pow(t,2)) )/ pow(t,2);

    double residue,r_min;
    bool loop;
    loop = true;

    while (loop == true){
        r_min = 0;

        for(int i =0; i< total_nodes;i++){
            if( ! bcs.get_bc(i)){
                residue = w*( 0.25* ( streamfunction[mesh.get_w_node(i)] + streamfunction[mesh.get_n_node(i)]+
                            streamfunction[mesh.get_s_node(i)] + streamfunction[mesh.get_e_node(i)]
                            + vorticity[i]) -streamfunction[i]);

                    r_min = r_min + fabs(residue);
                    streamfunction[i] = streamfunction[i] + residue;
            }else{

                streamfunction[i] = 0;
            }
        }
        r_min = r_min/total_nodes;
        std::cout << "R-min:" << r_min << std::endl;

     if (r_min < 10e-9){
        loop = false;
     }


    }


}
