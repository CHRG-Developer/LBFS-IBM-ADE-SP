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

using namespace std;
Solver::Solver()
{
    //ctor
}

Solver::~Solver()
{
    //dtor
}


void Solver::cell_interface_initialiser( double &rho_interface,vector_var &rho_u_interface,
                                        flux_var &x_flux,flux_var &y_flux ){
    // initialise variables
     // add in reset function
    rho_interface = 0;

    rho_u_interface.x =0;
    rho_u_interface.y = 0;
    rho_u_interface.z = 0;

    x_flux.P = 0;
    x_flux.momentum_x =0;
    x_flux.momentum_y =0;
    x_flux.momentum_z =0;

    y_flux.P = 0;
    y_flux.momentum_x =0;
    y_flux.momentum_y =0;
    y_flux.momentum_z =0;

}


double Solver::feq_calc_incomp(double weight, vector_var e_alpha, vector_var u_lattice, double u_magnitude,
                        double cs, double rho_lattice, double rho_0, int k){
    double feq;


    feq = e_alpha.Dot_Product(u_lattice) *3.0 ;
    feq = feq + ( pow(e_alpha.Dot_Product(u_lattice),2)  - pow((u_magnitude* cs),2) )
    *4.5;
    feq= feq *weight *rho_0 ;
     feq = feq + weight *rho_lattice ;

    return feq;

}


double Solver::feq_calc(double weight, vector_var e_alpha, vector_var u_lattice, double u_magnitude,
                        double cs, double rho_lattice){
    double feq;


    feq = 1.0  ;
    feq = feq
        + e_alpha.Dot_Product(u_lattice) *3.0 ;
    feq = feq + ( pow(e_alpha.Dot_Product(u_lattice),2)  - pow((u_magnitude* cs),2) )
    *4.5;
    feq= feq *weight *rho_lattice ;

    return feq;

}


void Solver::General_Purpose_Solver_mk_i( unstructured_mesh &Mesh , Solution &soln, Boundary_Conditions &bcs,
                                   external_forces &source,global_variables &globals, domain_geometry &domain,
                                   initial_conditions &init_conds, unstructured_bcs &quad_bcs_orig, int mg,
                                   Solution &residual, int fmg, post_processing &pp)
{

    ///Declarations
    RungeKutta rk4;
    Solution temp_soln(Mesh.get_total_cells()); // intermediate solution for RK
    Solution soln_t0(Mesh.get_total_cells()); // solution at t0 in RK cycle
    Solution soln_t1(Mesh.get_total_cells());
    Solution residual_worker(Mesh.get_total_cells()); // stores residuals
    Solution vortex_error(Mesh.get_total_cells());
    Solution real_error (Mesh.get_total_cells());
	Solution wall_shear_stress(Mesh.get_n_wall_cells());
    gradients grads (Mesh.get_total_cells());
	Solution cfl_areas(Mesh.get_total_cells());


    flux_var RK;

    double delta_t = globals.time_marching_step;

	double *delta_t_local;
	int *delta_t_frequency;
	double *local_time;
	bool *calc_cell;
	double *local_viscosity;
	double *local_fneq;

	delta_t_local = new double[Mesh.get_n_cells()];
	if (delta_t_local == NULL) exit(1);
	delta_t_frequency = new int[Mesh.get_n_cells()];
	if (delta_t_frequency == NULL) exit(1);
	local_time = new double[Mesh.get_n_cells()];
	if (local_time == NULL) exit(1);
	calc_cell = new bool[Mesh.get_n_cells()];
	if (calc_cell == NULL) exit(1);
	local_viscosity = new double[Mesh.get_total_cells()];
	if (local_viscosity == NULL) exit(1);
	local_fneq = new double[Mesh.get_total_cells()];
	if (local_fneq == NULL) exit(1);


	std::fill_n(calc_cell, Mesh.get_n_cells() , true);


    double local_tolerance;
    double rho_interface;
    double interface_area;
    double feq_lattice [15];
    double u_lattice[15], v_lattice[15], w_lattice[15], rho_lattice[15];

    double lattice_weight [15];
    double time;

    double f1,f2,f3,f4;
	double output_residual_threshold = 0;
    double visc;


    double angular_freq, wom_cos,force;

    std::ofstream error_output , vortex_output , max_u, debug_log;
    std::string output_dir,decay_dir,max_u_dir;
    output_dir = globals.output_file +"/error.txt";
    vector_var cell_1, cell_2, interface_node, lattice_node, delta_u, delta_v ,delta_w,delta_rho;
    vector_var relative_interface;
    vector_var  vel_lattice,  rho_u_interface , u_interface;
    vector_var delta_u1, delta_v1 ,delta_w1,delta_rho1;
    vector_var cell_normal;
    vector_var flux_e_alpha [9];
    vector_var u,v,w,rho;
    std::vector<vector_var> e_alpha;
    std::vector<int> cell_nodes;

    // vector_var flux_e_alpha;
    residuals convergence_residual;
    flux_var x_flux , y_flux,z_flux;
    flux_var cell_flux ;

    flux_var debug [4] ,debug_flux[4],arti_debug [4];
    flux_var dbug [4];
    flux_var int_debug[4];

    bc_var bc;

    int neighbour;
    int timesteps;
    int i;


    //calculate timesteps
    int center_node;
	int wall;

    center_node = Mesh.get_centre_node();

    tecplot_output tecplot;

    ///Initialisations

    dt = domain.dt; // timestepping for streaming // non-dim equals 1
    c = 1; // assume lattice spacing is equal to streaming timestep
    cs = c/sqrt(3);
    visc = (globals.tau -0.5)/3 * domain.dt;

    soln_t1.clone(soln);
    local_tolerance = globals.tolerance;
     delta_t =1 ;
    timesteps = ceil( globals.simulation_length);
    output_dir = globals.output_file +"/error.txt";
    decay_dir = globals.output_file +"/vortex_error.txt";
    max_u_dir = globals.output_file +"/max_u.txt";
   // error_output.open("/home/brendan/Dropbox/PhD/Test Cases/Couette Flow/error.txt", ios::out);
    error_output.open(output_dir.c_str(), ios::out);
    output_dir = globals.output_file +"/residual_log.txt";
    debug_log.open(output_dir.c_str(), ios::out);
    vortex_output.open(decay_dir.c_str(), ios::out);
    max_u.open(max_u_dir.c_str(), ios::out);


    populate_e_alpha(e_alpha,lattice_weight,c,globals.PI,15);
    time =0;
    angular_freq = visc* pow(globals.womersley_no,2) / pow(Mesh.get_Y()/2,2);
    force = -init_conds.pressure_gradient ;

    // residual_factor = delta_t/ Mesh.get_s_area(0);

    // taylor vortex memory
    double td;

	//drag co-efficients
	double drag_t1, drag_t0;
    time = 0;

    neighbour =0;

    td = 100000000000000000;

    grads.pre_fill_LHS_and_RHS_matrix(bcs,Mesh,domain,soln,globals);

    debug_log << "t,rk,i,res_rho,res_u,res_v,res_w,x,y,z, dt,visc,rho,u,v,ux,uy,uz,vx,vy,vz" << endl ;
//


	std::clock_t time1, time2, time3, time4;
	double duration1, duration2, duration3;




	populate_cfl_areas(cfl_areas,Mesh);



    // loop in time
    for (int t= 0; t < timesteps; t++){
        // soln is the solution at the start of every
		// RK step.(rk = n) Temp_soln holds the values at end of
		// step.(rk = n+1)
        soln_t0.clone(soln_t1);    // soln_t0 holds macro variable solution at start of time step
        temp_soln.clone(soln);  //temp holds rho,u,v                      // t= 0, rk = 0

        convergence_residual.reset();

        //womersley flow peculiarities
        if (globals.testcase == 4){
            wom_cos = cos(angular_freq * t * delta_t) ;
            force = -init_conds.pressure_gradient * wom_cos;
        }

        //find_real_time(delta_t_local, local_time, calc_face,Mesh,calc_cell);
        //local timestepping calculation
        get_cfl(delta_t,temp_soln,Mesh,globals,delta_t_local, delta_t_frequency,cfl_areas);



        for( int rk = 0; rk < rk4.timesteps; rk++){
			//time1 = clock();
			drag_t1 = 0.0;
            //temp_soln.Initialise();

            //update temp_soln boundary conditions
             temp_soln.update_unstructured_bcs(bcs,Mesh,domain,t);

             residual_worker.Initialise(); //set to zeros
                         //get gradients for

			 //time2 = clock();
            grads.Get_LS_Gradients(bcs,Mesh,domain,temp_soln,globals);

			//time3 = clock();

			//std::cout << "CPU Cycles Gradients:" << double(time3 - time2) << std::endl;
			wall = 0;
             // loop through each cell and exclude the ghost cells
             //using n_cells here rather than total_cells
			for (int face=0 ; face < Mesh.get_n_faces() ; face ++) {

                i = Mesh.get_mesh_owner(face);

                // volume initialisers
                interface_area = 0.0;

                //get current cell centre
                cell_1.x = Mesh.get_centroid_x(i);
                cell_1.y = Mesh.get_centroid_y(i);
                cell_1.z = Mesh.get_centroid_z(i);
                // add in reset function
                cell_flux.zero();

                //maybe take out m1 and m2 and calc outside loop

                cell_interface_variables( face, i,interface_node, neighbour,
                                         interface_area,cell_normal, bcs, bc, Mesh,
                                         cell_2,cell_1);



				if (face > Mesh.get_n_neighbours() && bcs.get_name(face - Mesh.get_n_neighbours()) == "empty") {
					interface_area = 0.0;
				}
              else{

					cell_interface_initialiser( rho_interface, rho_u_interface, x_flux,y_flux);

				   // dt for the cell interface
					dt =  Mesh.get_delta_t_face(face);

					//populate macro variables
					populate_lattice_macros(u_lattice,v_lattice,w_lattice,rho_lattice,cell_1,cell_2,
							interface_node,i,neighbour,grads,temp_soln,Mesh,bcs,cell_normal);

					//get initial feqs
					for(int k = 0; k< 15; k ++){
						populate_feq(u_lattice,v_lattice,w_lattice,rho_lattice,lattice_weight,
							feq_lattice,k,globals);
					}

					// get macroscopic values at cell interface
					 rho_interface = feq_lattice[0] + feq_lattice[1] + feq_lattice[2]+ feq_lattice[3]
					 +feq_lattice[4] + feq_lattice[5]+ feq_lattice[6]+feq_lattice[7]+ feq_lattice[8]
					 +feq_lattice[9] + feq_lattice[10]+ feq_lattice[11]+feq_lattice[12]+ feq_lattice[13]
					 + feq_lattice[14];

					u_interface.x = 1/rho_interface * ( feq_lattice[1]- feq_lattice[2]
					+feq_lattice[7] - feq_lattice[8]
					 +feq_lattice[9] - feq_lattice[10]+ feq_lattice[11] -feq_lattice[12]- feq_lattice[13]
					 + feq_lattice[14]);

					u_interface.y = 1/rho_interface * ( feq_lattice[3]- feq_lattice[4]
					+feq_lattice[7] - feq_lattice[8]
					 +feq_lattice[9] - feq_lattice[10]- feq_lattice[11] +feq_lattice[12]+ feq_lattice[13]
					 - feq_lattice[14]);

					u_interface.z = 1/rho_interface * ( feq_lattice[5]- feq_lattice[6]
					+feq_lattice[7] - feq_lattice[8]
					 - feq_lattice[9] + feq_lattice[10] + feq_lattice[11] -feq_lattice[12] + feq_lattice[13]
					 - feq_lattice[14]);

					calculate_flux_at_interface(u_interface,dt,globals,local_viscosity,rho_interface,lattice_weight
					,cell_flux,i,feq_lattice,cell_normal,interface_area,local_fneq);

					// add density flux to current cell and neighbouring cell
					residual_worker.add_rho(i,-cell_flux.P /Mesh.get_cell_volume(i));
					residual_worker.add_rho(neighbour, +cell_flux.P/Mesh.get_cell_volume(neighbour));

					// add x momentum
					residual_worker.add_u(i,-cell_flux.momentum_x/Mesh.get_cell_volume(i));
					residual_worker.add_u(neighbour, +cell_flux.momentum_x/Mesh.get_cell_volume(neighbour));

					// add y momentum
					residual_worker.add_v(i,-cell_flux.momentum_y/Mesh.get_cell_volume(i));
					residual_worker.add_v(neighbour, cell_flux.momentum_y/Mesh.get_cell_volume(neighbour));

					  // add z momentum
	                residual_worker.add_w(i,-cell_flux.momentum_z/Mesh.get_cell_volume(i));
	                residual_worker.add_w(neighbour, cell_flux.momentum_z/Mesh.get_cell_volume(neighbour));


					//get values for wall shear stress and drag
					if (face > Mesh.get_n_neighbours() && bcs.get_name(face - Mesh.get_n_neighbours()) == "wall") {
						wall_shear_stress.set_rho(wall, rho_interface);
						wall_shear_stress.set_u(wall, cell_flux.momentum_x);
						drag_t1 = drag_t1 + cell_flux.momentum_x;
						wall_shear_stress.set_v(wall, cell_flux.momentum_y);
						wall_shear_stress.set_w(wall, cell_flux.momentum_z);
						wall = wall + 1;
					}

                }

        }


            // for( int i=0; i < Mesh.get_total_cells(); i++){
            //
            //     debug_log << t << ", "  << rk << ", " << i << ", " << residual_worker.get_rho(i) << ", " <<
            //     residual_worker.get_u(i) << ", " << residual_worker.get_v(i) << ", " << residual_worker.get_w(i)
            //     << ", " <<
            //     Mesh.get_centroid_x(i) << " , " << Mesh.get_centroid_y(i) << "," <<  Mesh.get_centroid_z(i) << "," <<
            //      delta_t_local[i]  << " , " << local_viscosity[i] << "," <<
            //     temp_soln.get_rho(i)<< "," << temp_soln.get_u(i) << " , " << temp_soln.get_v(i)<< " , " <<
            //     grads.get_u(i).x << " , " << grads.get_u(i).y << " , " <<  grads.get_u(i).z << " , " <<
            //     grads.get_v(i).x << " , " << grads.get_v(i).y << " , " <<  grads.get_v(i).z << " , " <<
            //     grads.get_w(i).x << " , " << grads.get_w(i).y << "," <<  grads.get_w(i).z
            //
            //     << endl;
            //
            // }



          //  residual_worker.remove_double_errors();

            //update RK values
            for( int i=0; i < Mesh.get_n_cells(); i++){
                if( calc_cell[i]){

                    // update intermediate macroscopic variables for next Runge Kutta Time Step
                    f1 = soln_t0.get_rho(i) + residual_worker.get_rho(i)*delta_t_local[i] *rk4.alpha[rk];
                    f2 = soln_t0.get_u(i) + (residual_worker.get_u(i)+force) *delta_t_local[i]*rk4.alpha[rk];
                    f3 = soln_t0.get_v(i) + residual_worker.get_v(i) *delta_t_local[i]*rk4.alpha[rk];
                     f4 = soln_t0.get_w(i) + residual_worker.get_w(i) *delta_t_local[i]*rk4.alpha[rk];

                      // change momentum to velocity
                    f2 = f2/f1;
                    f3 =f3/f1;
                    f4=f4/f1;

                     temp_soln.update(f1,f2,f3,f4, i);
                    //temp_soln.update(1.0,f2,0.0,0.0, i);

                    //add contributions to
                    soln_t1.add_rho(i, delta_t_local[i]* rk4.beta[rk] * residual_worker.get_rho(i));
                    soln_t1.add_u(i, delta_t_local[i]* rk4.beta[rk] * (residual_worker.get_u(i)+force));
                    soln_t1.add_v(i, delta_t_local[i]* rk4.beta[rk] * residual_worker.get_v(i));
                    soln_t1.add_w(i, delta_t_local[i]* rk4.beta[rk] * residual_worker.get_w(i));


                    f1 = soln_t1.get_rho(i);
                    f2 = soln_t1.get_u(i)/soln_t1.get_rho(i);
                    f3 = soln_t1.get_v(i)/soln_t1.get_rho(i);
                    f4= soln_t1.get_w(i)/soln_t1.get_rho(i);

                   soln.update(f1,f2,f3,f4, i);
                   //soln.update(1.0,f2,0.0,0.0, i);
                }

            }

			//time4 = clock();
			//std::cout << "CPU Cycles Full RK cycle:" << double(time4- time1) << std::endl;
        }

        for( int i = 0; i < Mesh.get_n_cells(); i++){
            if(calc_cell[i]){
                    convergence_residual.add_l2_norm_residuals(residual_worker,i);
                       //error checking
                    if (std::isnan(temp_soln.get_rho(i)) || std::isnan(temp_soln.get_u(i))) {
                                    if( mg == 0){
                                        error_output.close();
                                    }

                                    cout << "nan failure" <<endl;
                                    cout << t <<endl;
                                    cout << i <<endl;

									delete[] delta_t_local;
									delta_t_local = NULL;
									delete[] delta_t_frequency;
									delta_t_frequency = NULL;
									delete[] local_time;
									local_time = NULL;
									delete[] calc_cell;
									calc_cell = NULL;
									delete[] local_viscosity;
									local_viscosity = NULL;
									delete[] local_fneq;
									local_fneq = NULL;

                                    return;
                            }
                    if (temp_soln.get_rho(i)/init_conds.average_rho > 1000.0){
                        cout << "rho failure" <<endl;
                        cout << t <<endl;
                        cout << i <<endl;

                        tecplot.tecplot_output_unstructured_soln(globals,Mesh,soln,bcs,time,pp, residual_worker,delta_t_local, local_fneq);

						delete[] delta_t_local;
						delta_t_local = NULL;
						delete[] delta_t_frequency;
						delta_t_frequency = NULL;
						delete[] local_time;
						local_time = NULL;
						delete[] calc_cell;
						calc_cell = NULL;
						delete[] local_viscosity;
						local_viscosity = NULL;
						delete[] local_fneq;
						local_fneq = NULL;

                        return;
                    }
            }
        }

        //convergence_residual.ansys_5_iter_rms(t);
        convergence_residual.l2_norm_rms_moukallad(globals);
        time = t*delta_t;

        if( mg == 0 && t%globals.output_step == 1){

            error_output << t << ", "  << convergence_residual.max_error()   << ", " <<
            convergence_residual.rho_rms << ", " << convergence_residual.u_rms << ", " <<
            convergence_residual.v_rms << ", " <<
				convergence_residual.w_rms << " , FMG cycle: " << fmg << endl;
            cout << "time t=" << time  << " error e =" << convergence_residual.max_error()
            <<  " delta_t:" << delta_t <<std::endl;
            max_u << t << "," << soln.get_u(center_node) << "," << force << endl;
			cout << "drag: " << drag_t1 << endl;

			//only output at decreasing order of magnitudes - save space on hard drive
			if (convergence_residual.max_error() < pow(10, output_residual_threshold)) {
				tecplot.tecplot_output_unstructured_soln(globals, Mesh, soln, bcs, time, pp, residual_worker, delta_t_local, local_fneq);
				output_residual_threshold = output_residual_threshold - 1;
				soln.output(globals.output_file, globals, domain);
			}
            //soln.output_centrelines(globals.output_file,globals,Mesh,time);

        }

        if ( convergence_residual.max_error() < local_tolerance || time > td){
            if( mg == 0){

                cout << "convergence" <<endl;
                cout << "time t=" << time  << " error e =" << convergence_residual.max_error()
                    <<  " delta_t:" << delta_t <<std::endl;
                error_output.close();
                debug_log.close();
                vortex_output.close();
                max_u.close();

                // vortex calcs
                temp_soln.update_unstructured_bcs(bcs,Mesh,domain,t);
                grads.Get_LS_Gradients(bcs,Mesh,domain,temp_soln,globals);
                pp.cylinder_post_processing(Mesh,globals,grads,bcs,temp_soln,domain,wall_shear_stress);
               // pp.calc_vorticity(x_gradients,y_gradients);
                 //pp.calc_streamfunction(Mesh,globals,bcs);
                 tecplot.tecplot_output_unstructured_soln(globals,Mesh,soln,bcs,time,pp,residual_worker, delta_t_local, local_fneq);
                //soln.output_centrelines(globals.output_file,globals,Mesh,time);
            }

			delete[] delta_t_local;
			delta_t_local = NULL;
			delete[] delta_t_frequency;
			delta_t_frequency = NULL;
			delete[] local_time;
			local_time = NULL;
			delete[] calc_cell;
			calc_cell = NULL;
			delete[] local_viscosity;
			local_viscosity = NULL;
			delete[] local_fneq;
			local_fneq = NULL;

            return ;
        }


    }

//    pp.calc_vorticity(x_gradients,y_gradients);
    //pp.calc_streamfunction(Mesh,globals,bcs);

    cout << "out of time" <<endl;
    error_output.close();
    vortex_output.close();
    debug_log.close();
    max_u.close();
	pp.cylinder_post_processing(Mesh, globals, grads, bcs, temp_soln, domain, temp_soln);
    tecplot.tecplot_output_unstructured_soln(globals,Mesh,soln,bcs,time,pp,residual_worker, delta_t_local,local_fneq);

	delete[] delta_t_local;
	delta_t_local = NULL;
	delete[] delta_t_frequency;
	delta_t_frequency = NULL;
	delete[] local_time;
	local_time = NULL;
	delete[] calc_cell;
	calc_cell = NULL;
	delete[] local_viscosity;
	local_viscosity = NULL;
	delete[] local_fneq;
	local_fneq = NULL;


}


void Solver::get_weighted_average( gradients &grads, int i, int neighbour, double m1, double m2,
   vector_var &u, vector_var &v, vector_var &w, vector_var &rho, unstructured_mesh &mesh)
{
    double a,b ,x,y,z;

    //check for boundary condition

    //use boundary cell gradients as these are at cell face
    if( neighbour > mesh.get_n_cells()){
         x =  grads.get_u(neighbour).x ;
        y = grads.get_u(neighbour).y ;
        z =  grads.get_u(neighbour).z ;
        u.set_equal(x,y,z);

        x = grads.get_v(neighbour).x ;
        y = grads.get_v(neighbour).y ;
        z = grads.get_v(neighbour).z ;
        v.set_equal(x,y,z);


        x = grads.get_w(neighbour).x ;
        y = grads.get_w(neighbour).y ;
        z = grads.get_w(neighbour).z ;
        w.set_equal(x,y,z);


        x =  grads.get_rho(neighbour).x ;
        y =  grads.get_rho(neighbour).y ;
        z =  grads.get_rho(neighbour).z ;
        rho.set_equal(x,y,z);


    }else{

        a = m1 +m2;
        b = m2/a;
        a = m1/a;


        x = grads.get_u(i).x * a + grads.get_u(neighbour).x *b;
        y = grads.get_u(i).y * a + grads.get_u(neighbour).y *b;
        z = grads.get_u(i).z * a + grads.get_u(neighbour).z *b;
        u.set_equal(x,y,z);

        x = grads.get_v(i).x * a + grads.get_v(neighbour).x *b;
        y = grads.get_v(i).y * a + grads.get_v(neighbour).y *b;
        z = grads.get_v(i).z * a + grads.get_v(neighbour).z *b;
        v.set_equal(x,y,z);


        x = grads.get_w(i).x * a + grads.get_w(neighbour).x *b;
        y = grads.get_w(i).y * a + grads.get_w(neighbour).y *b;
        z = grads.get_w(i).z * a + grads.get_w(neighbour).z *b;
        w.set_equal(x,y,z);


        x = grads.get_rho(i).x * a + grads.get_rho(neighbour).x *b;
        y = grads.get_rho(i).y * a + grads.get_rho(neighbour).y *b;
        z = grads.get_rho(i).z * a + grads.get_rho(neighbour).z *b;
        rho.set_equal(x,y,z);
        }


}

void Solver::calculate_flux_at_interface(vector_var u_interface, double dt, global_variables &globals,
    double local_viscosity[], double rho_interface, double lattice_weight[], flux_var &cell_flux,int i,
    double feq_lattice[], vector_var cell_normal,double interface_area, double local_fneq[]){

    double uu2, vv2, ww2, u2v2w2,uu,vv, ww, uv,uw,vw,fneq_tau;
    double feq_interface[15];
    flux_var x_flux, y_flux,z_flux;

        uu2 = u_interface.x * u_interface.x/ globals.pre_conditioned_gamma;
        vv2 = u_interface.y * u_interface.y/ globals.pre_conditioned_gamma;
        ww2 = u_interface.z *u_interface.z/ globals.pre_conditioned_gamma;

        u2v2w2 = (uu2 + vv2 + ww2) * 1.5;

        uu =u_interface.x;
        vv = u_interface.y;
        ww = u_interface.z;

        uv = uu*vv*9.0/ globals.pre_conditioned_gamma;
        uw = uu*ww*9.0/ globals.pre_conditioned_gamma;
        vw = vv*ww *9.0/ globals.pre_conditioned_gamma;


        fneq_tau = (globals.visc *3/dt / globals.pre_conditioned_gamma);
        local_viscosity[i] = fneq_tau/3*dt;
		local_fneq[i] = fneq_tau;


        feq_interface[1] = lattice_weight[1] * rho_interface*
              (1.0+3.0*uu+ 4.5*uu2 - u2v2w2);
        feq_interface[1] = feq_interface[1]
            - fneq_tau * (feq_interface[1] - feq_lattice[1]);

        feq_interface[2] = lattice_weight[2] * rho_interface*
            (1.0-3.0*uu + 4.5*uu2 - u2v2w2);
        feq_interface[2] = feq_interface[2]
            - fneq_tau * (feq_interface[2] - feq_lattice[2]);

        feq_interface[3] = lattice_weight[3] * rho_interface*
            (1.0 +3.0*vv + 4.5*vv2 - u2v2w2);
        feq_interface[3] = feq_interface[3]
            - fneq_tau * (feq_interface[3] - feq_lattice[3]);

        feq_interface[4] = lattice_weight[4] * rho_interface*
            (1.0 -3.0*vv + 4.5*vv2 - u2v2w2);
        feq_interface[4] = feq_interface[4]
            - fneq_tau * (feq_interface[4] - feq_lattice[4]);

        feq_interface[5] = lattice_weight[5] * rho_interface*
              (1.0 +3.0*ww + 4.5*ww2 - u2v2w2);
        feq_interface[5] = feq_interface[5]
            - fneq_tau * (feq_interface[5] - feq_lattice[5]);

        feq_interface[6] = lattice_weight[6] * rho_interface*
            (1.0 -3.0*ww + 4.5*ww2 - u2v2w2);
        feq_interface[6] = feq_interface[6]
            -fneq_tau * (feq_interface[6] - feq_lattice[6]);

        feq_interface[7] = lattice_weight[7] * rho_interface*
            (1.0 +3.0*uu +3.0*vv + 3.0*ww  +uv + uw + vw +
                4.5*uu2 + 4.5*vv2 + 4.5*ww2 -u2v2w2);
        feq_interface[7] = feq_interface[7]
            - fneq_tau * (feq_interface[7] - feq_lattice[7]);

        feq_interface[8] = lattice_weight[8] * rho_interface*
            (1.0 -3.0*uu -3.0*vv - 3.0*ww  +uv + uw + vw +
                4.5*uu2 + 4.5*vv2 + 4.5*ww2 -u2v2w2);
        feq_interface[8] = feq_interface[8]
            - fneq_tau * (feq_interface[8] - feq_lattice[8]);

         feq_interface[9] = lattice_weight[9] * rho_interface*
              (1.0 +3.0*uu +3.0*vv - 3.0*ww  +uv - uw - vw +
                4.5*uu2 + 4.5*vv2 + 4.5*ww2 -u2v2w2);
        feq_interface[9] = feq_interface[9]
            - fneq_tau * (feq_interface[9] - feq_lattice[9]);

         feq_interface[10] = lattice_weight[10] * rho_interface*
            (1.0 -3.0*uu -3.0*vv + 3.0*ww  +uv - uw - vw +
            4.5*uu2 + 4.5*vv2 + 4.5*ww2 -u2v2w2);
        feq_interface[10] = feq_interface[10]
            - fneq_tau * (feq_interface[10] - feq_lattice[10]);

         feq_interface[11] = lattice_weight[11] * rho_interface*
            (1.0 +3.0*uu -3.0*vv + 3.0*ww  -uv + uw - vw +
                4.5*uu2 + 4.5*vv2 + 4.5*ww2 -u2v2w2);
        feq_interface[11] = feq_interface[11]
            - fneq_tau * (feq_interface[11] - feq_lattice[11]);

         feq_interface[12] = lattice_weight[12] * rho_interface*
            (1.0 -3.0*uu +3.0*vv - 3.0*ww  -uv + uw - vw +
                4.5*uu2 + 4.5*vv2 + 4.5*ww2 -u2v2w2);
        feq_interface[12] = feq_interface[12]
            - fneq_tau * (feq_interface[12] - feq_lattice[12]);

         feq_interface[13] = lattice_weight[13] * rho_interface*
            (1.0 -3.0*uu +3.0*vv + 3.0*ww  -uv - uw + vw +
                4.5*uu2 + 4.5*vv2 + 4.5*ww2 -u2v2w2);
        feq_interface[13] = feq_interface[13]
            - fneq_tau * (feq_interface[13] - feq_lattice[13]);

         feq_interface[14] = lattice_weight[14] * rho_interface*
              (1.0 +3.0*uu -3.0*vv - 3.0*ww  -uv - uw + vw +
                4.5*uu2 + 4.5*vv2 + 4.5*ww2 -u2v2w2);
        feq_interface[14] = feq_interface[14]
            - fneq_tau * (feq_interface[14] - feq_lattice[14]);


        x_flux.P = ( feq_interface[1]- feq_interface[2]
                           +feq_interface[7] - feq_interface[8]
                            +feq_interface[9] - feq_interface[10]+ feq_interface[11]
                            -feq_interface[12]- feq_interface[13]
                            + feq_interface[14]);

        x_flux.momentum_x  =
            ( feq_interface[1]+ feq_interface[2]
                           +feq_interface[7] + feq_interface[8]
                            +feq_interface[9] + feq_interface[10]+ feq_interface[11]
                            +feq_interface[12]+ feq_interface[13]
                            + feq_interface[14]);

         x_flux.momentum_y  =
            feq_interface[7] + feq_interface[8]
                            +feq_interface[9] + feq_interface[10] - feq_interface[11] - feq_interface[12]
                            - feq_interface[13]
                            - feq_interface[14];

        x_flux.momentum_z  =
            feq_interface[7] + feq_interface[8]
                            -feq_interface[9] - feq_interface[10] + feq_interface[11] + feq_interface[12]
                            - feq_interface[13]
                            - feq_interface[14];



        y_flux.P = ( feq_interface[3]- feq_interface[4]
                           +feq_interface[7] - feq_interface[8]
                            +feq_interface[9] - feq_interface[10] - feq_interface[11]
                            + feq_interface[12] + feq_interface[13]
                            - feq_interface[14]);


        y_flux.momentum_x  =x_flux.momentum_y ;

        y_flux.momentum_y = ( feq_interface[3]+ feq_interface[4]
                           +feq_interface[7] + feq_interface[8]
                            +feq_interface[9] + feq_interface[10] + feq_interface[11]
                            + feq_interface[12] + feq_interface[13]
                            + feq_interface[14]);

        y_flux.momentum_z =  feq_interface[7] + feq_interface[8]
                            - feq_interface[9] - feq_interface[10] - feq_interface[11] - feq_interface[12]
                            + feq_interface[13]
                            + feq_interface[14];

         z_flux.P = ( feq_interface[5]- feq_interface[6]
                           +feq_interface[7] - feq_interface[8]
                            -feq_interface[9] + feq_interface[10] + feq_interface[11]
                            - feq_interface[12] + feq_interface[13]
                            - feq_interface[14]);
        z_flux.momentum_x  = x_flux.momentum_z;
        z_flux.momentum_y = y_flux.momentum_z;
        z_flux.momentum_z = ( feq_interface[5]+ feq_interface[6]
                           +feq_interface[7] + feq_interface[8]
                            +feq_interface[9] + feq_interface[10] + feq_interface[11]
                            + feq_interface[12] + feq_interface[13]
                            + feq_interface[14]);


        cell_flux.P =  (x_flux.P*cell_normal.x + y_flux.P * cell_normal.y +
                        z_flux.P *cell_normal.z)*interface_area ;
        cell_flux.momentum_x = (x_flux.momentum_x*cell_normal.x +
                                y_flux.momentum_x * cell_normal.y +
                                z_flux.momentum_x*cell_normal.z)*interface_area ;

        cell_flux.momentum_y = (x_flux.momentum_y*cell_normal.x +
                                y_flux.momentum_y * cell_normal.y +
                                z_flux.momentum_y*cell_normal.z)*interface_area ;


        cell_flux.momentum_z = (x_flux.momentum_z*cell_normal.x +
                                y_flux.momentum_z * cell_normal.y +
                                z_flux.momentum_z*cell_normal.z)*interface_area ;
}


void Solver::populate_lattice_macros(double u_lattice[],double v_lattice[],
            double w_lattice[],double rho_lattice[],vector_var cell_1, vector_var cell_2,
            vector_var interface_node, int i, int neighbour, gradients &grads, Solution &temp_soln,
           unstructured_mesh &mesh, Boundary_Conditions &bcs,vector_var &cell_normal){

            vector_var temp1, temp2,vel,upwind_temp;

            ///   case 0: // center node

            vel.set_equal(temp_soln.get_u(i), temp_soln.get_v(i), temp_soln.get_v(i));


            temp1 = cell_1;
            temp1.subtract(interface_node);
            temp2 = cell_2;
            temp2.subtract(interface_node);
            double rho_i, rho_nb, u_i,u_nb, v_i,v_nb, w_i,w_nb;

            int nb;
            nb = neighbour - mesh.get_n_cells();

            if( neighbour > mesh.get_n_cells()){
                // vboundary = v_i + grad_boundary * distance _ib

                if(bcs.get_rho_type(nb) == 1){
                    rho_lattice[0] = bcs.get_rho(nb);

                }else if(bcs.get_rho_type(nb) == 8){
                    rho_lattice[0] = temp_soln.get_rho(i) ;

                }else{
                    rho_lattice[0]  = temp_soln.get_rho(i) - temp1.Dot_Product(grads.get_rho(neighbour));
                }

                if(bcs.get_vel_type(nb) ==1){
                    u_lattice[0]  = bcs.get_u(nb);
                    v_lattice[0]  = bcs.get_v(nb);
                    w_lattice[0]  = bcs.get_w(nb);

                }else if(bcs.get_vel_type(nb) ==8){
                    u_lattice[0]  = temp_soln.get_u(i) ;
                    v_lattice[0]  = 0;
                    w_lattice[0]  = temp_soln.get_w(i)  ;

                }else{
                    u_lattice[0]  = temp_soln.get_u(i) - temp1.Dot_Product(grads.get_u(neighbour));
                    v_lattice[0]  = temp_soln.get_v(i) - temp1.Dot_Product(grads.get_v(neighbour));
                    w_lattice[0]  = temp_soln.get_w(i) - temp1.Dot_Product(grads.get_w(neighbour));
                }

                rho_i = rho_lattice[0] ;
                rho_nb = rho_lattice[0] ;

                u_i = u_lattice[0];
                u_nb = u_lattice[0];

                v_i = v_lattice[0];
                v_nb = v_lattice[0];

                w_i = w_lattice[0];
                w_nb = w_lattice[0];


            }else{

                rho_i = temp_soln.get_rho(i) - temp1.Dot_Product(grads.get_rho(i)) ;
                rho_nb  =  temp_soln.get_rho(neighbour) - temp2.Dot_Product(grads.get_rho(neighbour));
                rho_lattice[0]  = (rho_i + rho_nb)*0.5;

                u_i = temp_soln.get_u(i) - temp1.Dot_Product(grads.get_u(i)) ;
                u_nb  =  temp_soln.get_u(neighbour) - temp2.Dot_Product(grads.get_u(neighbour));
                u_lattice[0]  = (u_i + u_nb)*0.5;

                v_i = temp_soln.get_v(i) - temp1.Dot_Product(grads.get_v(i)) ;
                v_nb  =  temp_soln.get_v(neighbour) - temp2.Dot_Product(grads.get_v(neighbour));
                v_lattice[0]  = (v_i + v_nb)*0.5;

                w_i = temp_soln.get_w(i) - temp1.Dot_Product(grads.get_w(i)) ;
                w_nb  =  temp_soln.get_w(neighbour) - temp2.Dot_Product(grads.get_w(neighbour));
                w_lattice[0]  = (w_i + w_nb)*0.5;




            }


      ///  case 1:west_node

            if( -1* cell_normal.x > 0){
                rho_lattice[1] =  rho_nb - grads.get_rho(neighbour).x* dt;
                u_lattice[1] = u_nb- grads.get_u(neighbour).x* dt;
                v_lattice[1] = v_nb - grads.get_v(neighbour).x *dt;
                w_lattice[1] = w_nb - grads.get_w(neighbour).x *dt;

            }else{
                rho_lattice[1] =  rho_i - grads.get_rho(i).x* dt;
                u_lattice[1] = u_i- grads.get_u(i).x* dt;
                v_lattice[1] = v_i - grads.get_v(i).x *dt;
                w_lattice[1] = w_i - grads.get_w(i).x *dt;
            }

      ///  case 2: // east_node
             if( cell_normal.x > 0){
                rho_lattice[2] =  rho_nb + grads.get_rho(neighbour).x* dt;
                u_lattice[2] = u_nb+ grads.get_u(neighbour).x* dt;
                v_lattice[2] = v_nb + grads.get_v(neighbour).x *dt;
                w_lattice[2] = w_nb + grads.get_w(neighbour).x *dt;

            }else{
                rho_lattice[2] =  rho_i + grads.get_rho(i).x* dt;
                u_lattice[2] = u_i+ grads.get_u(i).x* dt;
                v_lattice[2] = v_i + grads.get_v(i).x *dt;
                w_lattice[2] = w_i + grads.get_w(i).x *dt;
            }

         ///   case 3: // bottom node

            if( -1* cell_normal.y > 0){
                rho_lattice[3] =  rho_nb - grads.get_rho(neighbour).y* dt;
                u_lattice[3] = u_nb- grads.get_u(neighbour).y* dt;
                v_lattice[3] = v_nb - grads.get_v(neighbour).y *dt;
                w_lattice[3] = w_nb - grads.get_w(neighbour).y *dt;

            }else{
                rho_lattice[3] =  rho_i - grads.get_rho(i).y* dt;
                u_lattice[3] = u_i- grads.get_u(i).y* dt;
                v_lattice[3] = v_i - grads.get_v(i).y *dt;
                w_lattice[3] = w_i - grads.get_w(i).y *dt;
            }


     ///   case 4: // top node

            if( cell_normal.y > 0){
                rho_lattice[4] =  rho_nb + grads.get_rho(neighbour).y* dt;
                u_lattice[4] = u_nb+ grads.get_u(neighbour).y* dt;
                v_lattice[4] = v_nb + grads.get_v(neighbour).y *dt;
                w_lattice[4] = w_nb + grads.get_w(neighbour).y *dt;

            }else{
                rho_lattice[4] =  rho_i + grads.get_rho(i).y* dt;
                u_lattice[4] = u_i+ grads.get_u(i).y* dt;
                v_lattice[4] = v_i + grads.get_v(i).y *dt;
                w_lattice[4] = w_i + grads.get_w(i).y *dt;
            }


    ///   case 5: // back node
            if( -1* cell_normal.z > 0){
                rho_lattice[5] =  rho_nb - grads.get_rho(neighbour).z* dt;
                u_lattice[5] = u_nb- grads.get_u(neighbour).z* dt;
                v_lattice[5] = v_nb - grads.get_v(neighbour).z *dt;
                w_lattice[5] = w_nb - grads.get_w(neighbour).z *dt;

            }else{
                rho_lattice[5] =  rho_i - grads.get_rho(i).z* dt;
                u_lattice[5] = u_i- grads.get_u(i).z* dt;
                v_lattice[5] = v_i - grads.get_v(i).z *dt;
                w_lattice[5] = w_i - grads.get_w(i).z *dt;
            }

    ///   case 6: // front node
            if( +1* cell_normal.z > 0){
                rho_lattice[6] =  rho_nb + grads.get_rho(neighbour).z* dt;
                u_lattice[6] = u_nb+ grads.get_u(neighbour).z* dt;
                v_lattice[6] = v_nb + grads.get_v(neighbour).z *dt;
                w_lattice[6] = w_nb + grads.get_w(neighbour).z *dt;

            }else{
                rho_lattice[6] =  rho_i + grads.get_rho(i).z* dt;
                u_lattice[6] = u_i + grads.get_u(i).z* dt;
                v_lattice[6] = v_i + grads.get_v(i).z *dt;
                w_lattice[6] = w_i + grads.get_w(i).z *dt;
            }


       /// case 7: back bottom west
            if( (-1* cell_normal.x + -1*cell_normal.y  + -1*cell_normal.z ) > 0){
                rho_lattice[7] =  rho_nb - grads.get_rho(neighbour).x* dt
                                         - grads.get_rho(neighbour).y* dt
                                         - grads.get_rho(neighbour).z* dt;
                u_lattice[7] =  u_nb - grads.get_u(neighbour).x* dt
                                         - grads.get_u(neighbour).y* dt
                                         - grads.get_u(neighbour).z* dt;
                v_lattice[7] =  v_nb - grads.get_v(neighbour).x* dt
                                         - grads.get_v(neighbour).y* dt
                                         - grads.get_v(neighbour).z* dt;
                w_lattice[7] =  w_nb - grads.get_w(neighbour).x* dt
                                         - grads.get_w(neighbour).y* dt
                                         - grads.get_w(neighbour).z* dt;

            }else{
                rho_lattice[7] =  rho_i - grads.get_rho(i).x* dt
                                         - grads.get_rho(i).y* dt
                                         - grads.get_rho(i).z* dt;
                u_lattice[7] =  u_i - grads.get_u(i).x* dt
                                         - grads.get_u(i).y* dt
                                         - grads.get_u(i).z* dt;
                v_lattice[7] =  v_i - grads.get_v(i).x* dt
                                         - grads.get_v(i).y* dt
                                         - grads.get_v(i).z* dt;
                w_lattice[7] =  w_i - grads.get_w(i).x* dt
                                         - grads.get_w(i).y* dt
                                         - grads.get_w(i).z* dt;
            }




       /// case 9: front bottom west
             if( (-1* cell_normal.x + -1*cell_normal.y  + 1*cell_normal.z ) > 0){
                rho_lattice[9] =  rho_nb - grads.get_rho(neighbour).x* dt
                                         - grads.get_rho(neighbour).y* dt
                                         + grads.get_rho(neighbour).z* dt;
                u_lattice[9] =  u_nb - grads.get_u(neighbour).x* dt
                                         - grads.get_u(neighbour).y* dt
                                         + grads.get_u(neighbour).z* dt;
                v_lattice[9] =  v_nb - grads.get_v(neighbour).x* dt
                                         - grads.get_v(neighbour).y* dt
                                         + grads.get_v(neighbour).z* dt;
                w_lattice[9] =  w_nb - grads.get_w(neighbour).x* dt
                                         - grads.get_w(neighbour).y* dt
                                         + grads.get_w(neighbour).z* dt;

            }else{
                rho_lattice[9] =  rho_i - grads.get_rho(i).x* dt
                                         - grads.get_rho(i).y* dt
                                         + grads.get_rho(i).z* dt;
                u_lattice[9] =  u_i - grads.get_u(i).x* dt
                                         - grads.get_u(i).y* dt
                                         + grads.get_u(i).z* dt;
                v_lattice[9] =  v_i - grads.get_v(i).x* dt
                                         - grads.get_v(i).y* dt
                                         + grads.get_v(i).z* dt;
                w_lattice[9] =  w_i - grads.get_w(i).x* dt
                                         - grads.get_w(i).y* dt
                                         + grads.get_w(i).z* dt;
            }



      ///  case 11: back top west
          if( (-1* cell_normal.x + cell_normal.y  + -1*cell_normal.z ) > 0){
                rho_lattice[11] =  rho_nb - grads.get_rho(neighbour).x* dt
                                         + grads.get_rho(neighbour).y* dt
                                         - grads.get_rho(neighbour).z* dt;
                u_lattice[11] =  u_nb - grads.get_u(neighbour).x* dt
                                         + grads.get_u(neighbour).y* dt
                                         - grads.get_u(neighbour).z* dt;
                v_lattice[11] =  v_nb - grads.get_v(neighbour).x* dt
                                         + grads.get_v(neighbour).y* dt
                                         - grads.get_v(neighbour).z* dt;
                w_lattice[11] =  w_nb - grads.get_w(neighbour).x* dt
                                         + grads.get_w(neighbour).y* dt
                                         - grads.get_w(neighbour).z* dt;

            }else{
                rho_lattice[11] =  rho_i - grads.get_rho(i).x* dt
                                         + grads.get_rho(i).y* dt
                                         - grads.get_rho(i).z* dt;
                u_lattice[11] =  u_i - grads.get_u(i).x* dt
                                         + grads.get_u(i).y* dt
                                         - grads.get_u(i).z* dt;
                v_lattice[11] =  v_i - grads.get_v(i).x* dt
                                         + grads.get_v(i).y* dt
                                         - grads.get_v(i).z* dt;
                w_lattice[11] =  w_i - grads.get_w(i).x* dt
                                         + grads.get_w(i).y* dt
                                         - grads.get_w(i).z* dt;
            }



      /// case 14: front top west
           if( (-1* cell_normal.x + cell_normal.y  + cell_normal.z ) > 0){
                rho_lattice[14] =  rho_nb - grads.get_rho(neighbour).x* dt
                                         + grads.get_rho(neighbour).y* dt
                                         + grads.get_rho(neighbour).z* dt;
                u_lattice[14] =  u_nb - grads.get_u(neighbour).x* dt
                                         + grads.get_u(neighbour).y* dt
                                         + grads.get_u(neighbour).z* dt;
                v_lattice[14] =  v_nb - grads.get_v(neighbour).x* dt
                                         + grads.get_v(neighbour).y* dt
                                         + grads.get_v(neighbour).z* dt;
                w_lattice[14] =  w_nb - grads.get_w(neighbour).x* dt
                                         + grads.get_w(neighbour).y* dt
                                         + grads.get_w(neighbour).z* dt;

            }else{
                rho_lattice[14] =  rho_i - grads.get_rho(i).x* dt
                                         + grads.get_rho(i).y* dt
                                         + grads.get_rho(i).z* dt;
                u_lattice[14] =  u_i - grads.get_u(i).x* dt
                                         + grads.get_u(i).y* dt
                                         + grads.get_u(i).z* dt;
                v_lattice[14] =  v_i - grads.get_v(i).x* dt
                                         + grads.get_v(i).y* dt
                                         + grads.get_v(i).z* dt;
                w_lattice[14] =  w_i - grads.get_w(i).x* dt
                                         + grads.get_w(i).y* dt
                                         + grads.get_w(i).z* dt;
            }



       /// case 8: front top east
         if( (1* cell_normal.x + 1*cell_normal.y  + 1*cell_normal.z ) > 0){
                rho_lattice[8] =  rho_nb + grads.get_rho(neighbour).x* dt
                                         + grads.get_rho(neighbour).y* dt
                                         + grads.get_rho(neighbour).z* dt;
                u_lattice[8] =  u_nb + grads.get_u(neighbour).x* dt
                                         + grads.get_u(neighbour).y* dt
                                         + grads.get_u(neighbour).z* dt;
                v_lattice[8] =  v_nb + grads.get_v(neighbour).x* dt
                                         + grads.get_v(neighbour).y* dt
                                         + grads.get_v(neighbour).z* dt;
                w_lattice[8] =  w_nb + grads.get_w(neighbour).x* dt
                                         + grads.get_w(neighbour).y* dt
                                         + grads.get_w(neighbour).z* dt;

            }else{
                rho_lattice[8] =  rho_i + grads.get_rho(i).x* dt
                                         + grads.get_rho(i).y* dt
                                         + grads.get_rho(i).z* dt;
                u_lattice[8] =  u_i + grads.get_u(i).x* dt
                                         + grads.get_u(i).y* dt
                                         + grads.get_u(i).z* dt;
                v_lattice[8] =  v_i + grads.get_v(i).x* dt
                                         + grads.get_v(i).y* dt
                                         + grads.get_v(i).z* dt;
                w_lattice[8] =  w_i + grads.get_w(i).x* dt
                                         + grads.get_w(i).y* dt
                                         + grads.get_w(i).z* dt;
            }





         /// case 10 Back Top East
         if( (1* cell_normal.x + 1*cell_normal.y  + -1*cell_normal.z ) > 0){
                rho_lattice[10] =  rho_nb + grads.get_rho(neighbour).x* dt
                                         + grads.get_rho(neighbour).y* dt
                                         - grads.get_rho(neighbour).z* dt;
                u_lattice[10] =  u_nb + grads.get_u(neighbour).x* dt
                                         + grads.get_u(neighbour).y* dt
                                         - grads.get_u(neighbour).z* dt;
                v_lattice[10] =  v_nb + grads.get_v(neighbour).x* dt
                                         + grads.get_v(neighbour).y* dt
                                         - grads.get_v(neighbour).z* dt;
                w_lattice[10] =  w_nb + grads.get_w(neighbour).x* dt
                                         + grads.get_w(neighbour).y* dt
                                         - grads.get_w(neighbour).z* dt;

            }else{
                rho_lattice[10] =  rho_i + grads.get_rho(i).x* dt
                                         + grads.get_rho(i).y* dt
                                         - grads.get_rho(i).z* dt;
                u_lattice[10] =  u_i + grads.get_u(i).x* dt
                                         + grads.get_u(i).y* dt
                                         - grads.get_u(i).z* dt;
                v_lattice[10] =  v_i + grads.get_v(i).x* dt
                                         + grads.get_v(i).y* dt
                                         - grads.get_v(i).z* dt;
                w_lattice[10] =  w_i + grads.get_w(i).x* dt
                                         + grads.get_w(i).y* dt
                                         - grads.get_w(i).z* dt;
            }



         /// case 12 Front Bottom East
         if( (1* cell_normal.x + -1*cell_normal.y  + 1*cell_normal.z ) > 0){
                rho_lattice[12] =  rho_nb + grads.get_rho(neighbour).x* dt
                                         - grads.get_rho(neighbour).y* dt
                                         + grads.get_rho(neighbour).z* dt;
                u_lattice[12] =  u_nb + grads.get_u(neighbour).x* dt
                                         - grads.get_u(neighbour).y* dt
                                         + grads.get_u(neighbour).z* dt;
                v_lattice[12] =  v_nb + grads.get_v(neighbour).x* dt
                                         - grads.get_v(neighbour).y* dt
                                         + grads.get_v(neighbour).z* dt;
                w_lattice[12] =  w_nb + grads.get_w(neighbour).x* dt
                                         - grads.get_w(neighbour).y* dt
                                         + grads.get_w(neighbour).z* dt;

            }else{
                rho_lattice[12] =  rho_i + grads.get_rho(i).x* dt
                                         - grads.get_rho(i).y* dt
                                         + grads.get_rho(i).z* dt;
                u_lattice[12] =  u_i + grads.get_u(i).x* dt
                                         - grads.get_u(i).y* dt
                                         + grads.get_u(i).z* dt;
                v_lattice[12] =  v_i + grads.get_v(i).x* dt
                                         - grads.get_v(i).y* dt
                                         + grads.get_v(i).z* dt;
                w_lattice[12] =  w_i + grads.get_w(i).x* dt
                                         - grads.get_w(i).y* dt
                                         + grads.get_w(i).z* dt;
            }


         /// case 13 Back Bottom East
            if( (1* cell_normal.x + -1*cell_normal.y  + -1*cell_normal.z ) > 0){
                rho_lattice[13] =  rho_nb + grads.get_rho(neighbour).x* dt
                                         - grads.get_rho(neighbour).y* dt
                                         - grads.get_rho(neighbour).z* dt;
                u_lattice[13] =  u_nb + grads.get_u(neighbour).x* dt
                                         - grads.get_u(neighbour).y* dt
                                         - grads.get_u(neighbour).z* dt;
                v_lattice[13] =  v_nb + grads.get_v(neighbour).x* dt
                                         - grads.get_v(neighbour).y* dt
                                         - grads.get_v(neighbour).z* dt;
                w_lattice[13] =  w_nb + grads.get_w(neighbour).x* dt
                                         - grads.get_w(neighbour).y* dt
                                         - grads.get_w(neighbour).z* dt;

            }else{
                rho_lattice[13] =  rho_i + grads.get_rho(i).x* dt
                                         - grads.get_rho(i).y* dt
                                         - grads.get_rho(i).z* dt;
                u_lattice[13] =  u_i + grads.get_u(i).x* dt
                                         - grads.get_u(i).y* dt
                                         - grads.get_u(i).z* dt;
                v_lattice[13] =  v_i + grads.get_v(i).x* dt
                                         - grads.get_v(i).y* dt
                                         - grads.get_v(i).z* dt;
                w_lattice[13] =  w_i + grads.get_w(i).x* dt
                                         - grads.get_w(i).y* dt
                                         - grads.get_w(i).z* dt;
            }


}



void Solver::populate_feq(double u_lattice[],double v_lattice[],
            double w_lattice[],double rho_lattice[],double lattice_weight[],
            double feq_lattice[],int k,global_variables &globals){

    ///d3q15 velocity set

    double uu2, vv2,u2v2w2,uv,uu,vv,ww2,uw,vw,ww;

    uu2 = u_lattice[k] * u_lattice[k] / globals.pre_conditioned_gamma;
    vv2 = v_lattice[k]* v_lattice[k] / globals.pre_conditioned_gamma;
    ww2 = w_lattice[k] *w_lattice[k] / globals.pre_conditioned_gamma;
    u2v2w2 = (uu2 + vv2 + ww2) * 1.5 ;

    uv = u_lattice[k]*v_lattice[k]*9.0 / globals.pre_conditioned_gamma;
    uw = u_lattice[k]*w_lattice[k]*9.0 / globals.pre_conditioned_gamma;
    vw = v_lattice[k]*w_lattice[k]*9.0 / globals.pre_conditioned_gamma;

    uu =u_lattice[k];
    vv = v_lattice[k];
    ww = w_lattice[k];

    switch(k) {

    case 0:
         feq_lattice[k] = lattice_weight[k] * rho_lattice[k]*(1.0

        -u2v2w2)  ;
        break;

    case 1:
        feq_lattice[k] = lattice_weight[k] * rho_lattice[k]*
        (1.0+3.0*uu+ 4.5*uu2 - u2v2w2);
        break;
    case 2:
        feq_lattice[k] = lattice_weight[k] * rho_lattice[k]*
        (1.0-3.0*uu + 4.5*uu2 - u2v2w2);
        break;
    case 3:
        feq_lattice[k] = lattice_weight[k] * rho_lattice[k]*
        (1.0 +3.0*vv + 4.5*vv2 - u2v2w2);
        break;
    case 4:
        feq_lattice[k] = lattice_weight[k] * rho_lattice[k]*
        (1.0 -3.0*vv + 4.5*vv2 - u2v2w2);
        break;
    case 5:
        feq_lattice[k] = lattice_weight[k] * rho_lattice[k]*
        (1.0 +3.0*ww + 4.5*ww2 - u2v2w2);
        break;
   case 6:
        feq_lattice[k] = lattice_weight[k] * rho_lattice[k]*
        (1.0 -3.0*ww + 4.5*ww2 - u2v2w2);
        break;
    case 7:
        feq_lattice[k] = lattice_weight[k] * rho_lattice[k]*
        (1.0 +3.0*uu +3.0*vv + 3.0*ww  +uv + uw + vw +
        4.5*uu2 + 4.5*vv2 + 4.5*ww2 -u2v2w2);
        break;
    case 8:
        feq_lattice[k] = lattice_weight[k] * rho_lattice[k]*
        (1.0 -3.0*uu -3.0*vv - 3.0*ww  +uv + uw + vw +
        4.5*uu2 + 4.5*vv2 + 4.5*ww2 -u2v2w2);
        break;
    case 9:
        feq_lattice[k] = lattice_weight[k] * rho_lattice[k]*
        (1.0 +3.0*uu +3.0*vv - 3.0*ww  +uv - uw - vw +
        4.5*uu2 + 4.5*vv2 + 4.5*ww2 -u2v2w2);
        break;
    case 10:
        feq_lattice[k] = lattice_weight[k] * rho_lattice[k]*
        (1.0 -3.0*uu -3.0*vv + 3.0*ww  +uv - uw - vw +
        4.5*uu2 + 4.5*vv2 + 4.5*ww2 -u2v2w2);
        break;
    case 11:
        feq_lattice[k] = lattice_weight[k] * rho_lattice[k]*
        (1.0 +3.0*uu -3.0*vv + 3.0*ww  -uv + uw - vw +
        4.5*uu2 + 4.5*vv2 + 4.5*ww2 -u2v2w2);
        break;
    case 12:
        feq_lattice[k] = lattice_weight[k] * rho_lattice[k]*
        (1.0 -3.0*uu +3.0*vv - 3.0*ww  -uv + uw - vw +
        4.5*uu2 + 4.5*vv2 + 4.5*ww2 -u2v2w2);
        break;
    case 13:
        feq_lattice[k] = lattice_weight[k] * rho_lattice[k]*
        (1.0 -3.0*uu +3.0*vv + 3.0*ww  -uv - uw + vw +
        4.5*uu2 + 4.5*vv2 + 4.5*ww2 -u2v2w2);
        break;
    case 14:
        feq_lattice[k] = lattice_weight[k] * rho_lattice[k]*
        (1.0 +3.0*uu -3.0*vv - 3.0*ww  -uv - uw + vw +
        4.5*uu2 + 4.5*vv2 + 4.5*ww2 -u2v2w2);
        break;
    }


}


vector_var Solver::get_e_alpha(int k, double &lattice_weight, double c, double PI ){

        vector_var temp;
        int x ,y,z;
        //get e_alpha again
        if (k >0 && k< 5){ //

            x = round(cos((k-1)*PI/2 ) * c );
            y = round(sin((k-1)*PI/2 )* c);
            z = 0; //update in 3D
            lattice_weight = 1.0/9.0;
        }else if( k >4){

            x = round(sqrt(2) * cos((k-5)*PI/2 + PI/4 ) * c );
            y = round(sqrt(2) * sin((k-5)*PI/2 + PI/4 ) * c);
            z = 0; //update in 3D
            lattice_weight = 1.0/36.0;

        }else{
            x = 0 ;
            y = 0;
            z = 0;
            lattice_weight = 4.0/9.0;
        }
        temp.x = x;
        temp.y = y;
        temp.z = z;


    return temp;
}

void Solver::populate_e_alpha(vector<vector_var> &e_alpha, double *lattice_weight, double c, double PI,int j ){

        vector_var temp;
        int x[15] = {0,1,-1,0,0,0,0,1,-1, 1,-1,1,-1,-1,1};
        int y[15] = { 0,0,0,1,-1,0,0,1,-1,1,-1,-1,1,1,-1};
        int z[15] = { 0,0,0,0,0,1,-1,1,-1,-1,1,1,-1,1,-1};
        //get e_alpha again

        for(int k =0; k<j; k++){
            if (k >0 && k< 7){ //

                lattice_weight[k] = 1.0/9.0;

            }else if( k >6){


                lattice_weight[k] = 1.0/72.0;

            }else{

                lattice_weight[k] = 2.0/9.0;
            }



            temp.x = x[k];
            temp.y = y[k];
            temp.z = z[k];

            e_alpha.push_back(temp);


        }



}

void Solver::get_cell_gradients(Mesh &Mesh, int i, int neighbour, int j, Solution &temp_soln,
                                vector_var &delta_rho, vector_var &delta_rho1,
                                vector_var &delta_u, vector_var &delta_u1,
                                vector_var &delta_v, vector_var &delta_v1,
                                Boundary_Conditions &bcs){


        int neighbour_1, neighbour_2;
        vector_var cell_1, cell_2;
        // is it N-S or E-W
        if( j == 2){


            neighbour_1 = Mesh.get_w_node(i);
            neighbour_2 = Mesh.get_e_node(i);

        }else{
            neighbour_1 = Mesh.get_s_node(i);
            neighbour_2 = Mesh.get_n_node(i);

        }

        // get neighbouring cells of cells
        Mesh.get_centroid(neighbour_1,cell_1);
        Mesh.get_centroid(neighbour_2,cell_2);

        delta_rho.Get_Gradient(temp_soln.get_rho(neighbour_1),temp_soln.get_rho(neighbour_2)
                               ,cell_1,cell_2);
        delta_u.Get_Gradient(temp_soln.get_u(neighbour_1),temp_soln.get_u(neighbour_2)
                               ,cell_1,cell_2);
        delta_v.Get_Gradient(temp_soln.get_v(neighbour_1),temp_soln.get_v(neighbour_2)
                               ,cell_1,cell_2);


        // get gradient of neighbouring cell
          if( j == 2){

            neighbour_1 = Mesh.get_w_node(neighbour);
            neighbour_2 = Mesh.get_e_node(neighbour);

        }else{
            neighbour_1 = Mesh.get_s_node(neighbour);
            neighbour_2 = Mesh.get_n_node(neighbour);

        }

        // get neighbouring cells of cells
        Mesh.get_centroid(neighbour_1,cell_1);
        Mesh.get_centroid(neighbour_2,cell_2);

        delta_rho1.Get_Gradient(temp_soln.get_rho(neighbour_1),temp_soln.get_rho(neighbour_2)
                               ,cell_1,cell_2);
        delta_u1.Get_Gradient(temp_soln.get_u(neighbour_1),temp_soln.get_u(neighbour_2)
                               ,cell_1,cell_2);
        delta_v1.Get_Gradient(temp_soln.get_v(neighbour_1),temp_soln.get_v(neighbour_2)
                               ,cell_1,cell_2);

        }

void Solver::cell_interface_variables( int j, int i, vector_var &interface_node, int &neighbour, double &interface_area,
                              vector_var &cell_normal, Boundary_Conditions &boundary_conditions,  bc_var &bc,
                              Mesh &Mesh, vector_var &cell_2) {

          switch(j) {

            case 0: // West
                interface_node.x = Mesh.get_west_x(i);
                interface_node.y = Mesh.get_west_y(i);
                interface_node.z= Mesh.get_west_z(i);
                neighbour = Mesh.get_w_node(i);
                interface_area = Mesh.get_w_area(i);
                cell_normal.x = Mesh.get_w_i(i);
                cell_normal.y = Mesh.get_w_j(i);
                cell_normal.z = Mesh.get_w_k(i);
                break;

            case 1: // South
                interface_node.x = Mesh.get_south_x(i);
                interface_node.y = Mesh.get_south_y(i);
                interface_node.z= Mesh.get_south_z(i);
                neighbour =Mesh.get_s_node(i);
                interface_area = Mesh.get_s_area(i);
                cell_normal.x = Mesh.get_s_i(i);
                cell_normal.y = Mesh.get_s_j(i);
                cell_normal.z = Mesh.get_s_k(i);

                break;
            case 2: // East
                interface_node.x = Mesh.get_east_x(i);
                interface_node.y = Mesh.get_east_y(i);
                interface_node.z= Mesh.get_east_z(i);
                interface_area = Mesh.get_e_area(i);
                neighbour =Mesh.get_e_node(i);
                cell_normal.x = Mesh.get_e_i(i);
                cell_normal.y = Mesh.get_e_j(i);
                cell_normal.z = Mesh.get_e_k(i);

                break;
            case 3: // North
                interface_node.x = Mesh.get_north_x(i);
                interface_node.y = Mesh.get_north_y(i);
                interface_node.z= Mesh.get_north_z(i);
                neighbour =Mesh.get_n_node(i);
                interface_area = Mesh.get_n_area(i);
                cell_normal.x = Mesh.get_n_i(i);
                cell_normal.y = Mesh.get_n_j(i);
                cell_normal.z = Mesh.get_n_k(i);

                break;
            case 4: // Front
                interface_node.x = Mesh.get_front_x(i);
                interface_node.y = Mesh.get_front_y(i);
                interface_node.z= Mesh.get_front_z(i);
                neighbour = Mesh.get_f_node(i);
                interface_area = Mesh.get_f_area(i);
                cell_normal.x = Mesh.get_f_i(i);
                cell_normal.y = Mesh.get_f_j(i);
                cell_normal.z = Mesh.get_f_k(i);

                break;
            case 5: // Back
                interface_node.x = Mesh.get_back_x(i);
                interface_node.y = Mesh.get_back_y(i);
                interface_node.z= Mesh.get_back_z(i);
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



      void Solver::cell_interface_variables( int face, int i, vector_var &interface_node, int &neighbour, double &interface_area,
                              vector_var &cell_normal, Boundary_Conditions &boundary_conditions,  bc_var &bc,
                              unstructured_mesh &Mesh, vector_var &cell_2, vector_var &cell_1 ) {



        interface_node.x = Mesh.get_face_x(face);
        interface_node.y = Mesh.get_face_y(face);
        interface_node.z= Mesh.get_face_z(face);

        neighbour = Mesh.get_mesh_neighbour(face);
        interface_area = Mesh.get_face_area(face);
        cell_normal.x = Mesh.get_face_i(face);
        cell_normal.y = Mesh.get_face_j(face);
        cell_normal.z = Mesh.get_face_k(face);


       cell_2.x = Mesh.get_centroid_x(neighbour);
       cell_2.y = Mesh.get_centroid_y((neighbour));
       cell_2.z = Mesh.get_centroid_z(neighbour);


      }



      void Solver::get_cell_nodes(std::vector<int> &cell_nodes, Boundary_Conditions &bcs,int neighbour,
                                Mesh &Mesh , int i,int j){

         //current cell
            cell_nodes.clear();
            if(bcs.get_bc(i) || bcs.get_bc(neighbour)){
                 cell_nodes.push_back(i);
                 cell_nodes.push_back(neighbour);

            }else if( j ==2){
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
            }else{
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
	  void Solver::populate_cfl_areas(Solution &cfl_areas, unstructured_mesh &Mesh ) {

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




      //get CFL numbers for inviscid and viscous matrices
      // see what time stepping results
void Solver::get_cfl( double &delta_t, Solution &soln
                ,unstructured_mesh &Mesh , global_variables &globals, double* delta_t_local, int* delta_t_frequency,Solution &cfl_areas){

    double cfl_inviscid[Mesh.get_total_cells()];
    double cfl_viscous[Mesh.get_total_cells()];
    double cfl_i_max,cfl_v_max;
    double t_v, t_i;
    double factor;
    cfl_i_max = 0;
    cfl_v_max = 0;
    int face,neighbour;

    double area_x_eigen,eigen,vel_mag,visc_eigen;
    factor = globals.time_marching_step;

    double visc_constant;
    visc_constant = 4;

    double min_delta_t,temp;

    double effective_speed_of_sound;
    //effective_speed_of_sound = 1/sqrt(3);
	  effective_speed_of_sound = globals.max_velocity* sqrt( 1 - globals.pre_conditioned_gamma + pow( globals.pre_conditioned_gamma /sqrt(3) / globals.max_velocity, 2));
    //loop through cells

    min_delta_t = 100000000000;

    for (int i =0; i< Mesh.get_n_cells(); i++){
		delta_t_frequency[i] = 1;

		// eigen values as per Zhaoli guo(2004) - preconditioning

		  //estimation of spectral radii s per Jiri Blasek: CFD Principles and Application Determination of Max time Step

        area_x_eigen = 0;
		area_x_eigen = (fabs(soln.get_u(i)) + effective_speed_of_sound)*cfl_areas.get_u(i)
			+ ( fabs(soln.get_v(i)) + effective_speed_of_sound)*cfl_areas.get_v(i)
			+ (fabs(soln.get_w(i)) + effective_speed_of_sound)*cfl_areas.get_w(i);

		area_x_eigen = area_x_eigen / globals.pre_conditioned_gamma;


		//reducing preconditioning increases viscous flux - increases eigenvalue
		visc_eigen = 2 * globals.visc / globals.pre_conditioned_gamma / soln.get_rho(i) / Mesh.get_cell_volume(i);
		visc_eigen = visc_eigen * (cfl_areas.get_u(i)*cfl_areas.get_u(i) + cfl_areas.get_v(i)*cfl_areas.get_v(i) + cfl_areas.get_w(i)* cfl_areas.get_w(i));

		area_x_eigen = area_x_eigen + 4 * visc_eigen;

        // use smallest time step allowed
        temp = factor*Mesh.get_cell_volume(i)/area_x_eigen;
		if (temp < 0) {
			min_delta_t = temp;


		}
		if (temp < min_delta_t) {
			min_delta_t = temp;
		}

        if (globals.time_stepping == "local" || globals.time_stepping == "talts" ){
            delta_t_local[i] = temp;

        }else{ //constant user defined time step
            delta_t_local[i] = factor;

		}
	}

	if( globals.time_stepping == "min"){
        std::fill_n(delta_t_local, Mesh.get_n_cells() , min_delta_t);
    }


	if (globals.time_stepping == "talts") {

		for (int i = 0; i < Mesh.get_n_cells(); i++) {
			delta_t_frequency[i] = pow(2, floor(log2(delta_t_local[i] / min_delta_t)) ) ;
			delta_t_local[i] = min_delta_t * delta_t_frequency[i] ;
		}

	}

    return;
}

void Solver::inverse_weighted_distance_interpolation(double &u, double &v, double &rho, Boundary_Conditions &bcs,
        Mesh &Mesh , domain_geometry &domain,Solution &soln,vector_var &interface_node,
         int k ,int i ,int neighbour,vector<vector_var> &e_alpha, int j,std::vector<int> &cell_nodes){

            // get interface node
            double w_u, w_v, w_rho, w_sum ,w;  // weighted macros

            // get 8 nodes'
            w_u = 0.0;
            w_v =0.0;
            w_rho =0.0;
            w_sum = 0.0;

            double r;
            r = 0.0;
            double dt;
            if (j == 2){
                dt = Mesh.get_delta_t_e(i);
            }else{
                dt = Mesh.get_delta_t_n(i);
            }

                      //get displacements
            vector_var node_displacement, target_node;

            // get target node
            target_node.x = interface_node.x -e_alpha[k].x * dt;
            target_node.y = interface_node.y -e_alpha[k].y * dt;
            target_node.z = interface_node.z -e_alpha[k].z * dt;

         for(auto &it : cell_nodes){
                node_displacement.x = Mesh.get_centroid_x(it) -target_node.x;
                node_displacement.y = Mesh.get_centroid_y(it) -target_node.y;
                node_displacement.z = Mesh.get_centroid_z(it) - target_node.z;

                r = node_displacement.Magnitude();

                //
                if (r < 10e-5){
                    u = soln.get_u(it);
                    w = soln.get_v(it);
                    rho = soln.get_rho(it);
                    return;

                }

                //get weight for this cc
                w = pow(1/r, 2.0);

                // sum weighted cc values
                w_u = w_u  + w* soln.get_u(it);
                w_v = w_v  + w* soln.get_v(it);
                w_rho = w_rho  + w* soln.get_rho(it);
                w_sum = w_sum + w;

            }

            // calc u v rho for target node
            u = w_u /w_sum;
            v = w_v/ w_sum;
            rho = w_rho/ w_sum;

        }

void Solver::find_real_time(double* delta_t_local, double* local_time, bool* calc_face,
                    unstructured_mesh &Mesh, bool* calc_cell){

    // for each cell check cell calc check if time is greater than neighbouring cells;
    int nb;

    for(int i =0; i<  Mesh.get_total_cells(); i ++){
        // first update actual time
        if(calc_cell[i]){
            local_time[i] = local_time[i] + delta_t_local[i];
        }
    }


    for( int i = 0; i < Mesh.get_total_cells() ; i++){


        calc_cell[i] = true;
        for(int j = 0; j < Mesh.gradient_cells[i].size(); j++){
           nb = Mesh.gradient_cells[i][j];
            if( local_time[i] > local_time[nb]){
                calc_cell[i] = false;
                j = Mesh.gradient_cells[i].size();
            }
        }
    }

    // then for each face calc if it should be calculated
    for( int k = 0; k < Mesh.get_n_faces(); k++){
        calc_face[k] = false;
        if( calc_cell[Mesh.get_mesh_owner(k)] || calc_cell[Mesh.get_mesh_neighbour(k)]){
            calc_face[k] = true;
        }
    }




}
