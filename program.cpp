#include "program.h"
#include "quad_bcs.h"
#include "global_variables.h"
#include "preprocessor.h"
#include <ctime>
#include "quad_bcs_plus.h"
#include "domain_geometry.h"
#include "Mesh.h"
#include "unstructured_mesh.h"
#include "Boundary_Conditions.h"
#include "external_forces.h"
#include "Solution.h"
#include "Solver.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "TECIO.h"
#include <sstream>
#include "tecplot_output.h"
#include <boost/filesystem.hpp>
#include "post_processing.h"
#include "unstructured_bcs.h"
#include "Solver_gpu.h"
#include "lagrangian_object.h"

namespace fs = boost::filesystem;

program::program()
{
    //ctor
}

program::~program()
{
    //dtor
}


void program::run(char* xml_input){


    quad_bcs_plus bcs;
    unstructured_bcs u_bcs;
    global_variables globals;
    domain_geometry domain;
    initial_conditions initial_conds;
	gpu_solver gpu;
     Solver solve;
    int mg =0; // first multigrid cycle
    int fmg =0; // first full multigrid cycle
    std::clock_t start = clock();
    double duration;
    duration = 0;
	std::vector< lagrangian_object> lagrangian_object_vec;

    preprocessor pre_processor;

    //postprocessor post_processor;
	std::cout << "Initialising Globals" << endl;

    pre_processor.initialise_program_variables(xml_input, globals, domain,initial_conds,bcs,u_bcs, lagrangian_object_vec);
	tecplot_output tp;

    copyfile(xml_input,globals.output_file);

    remove_existing_files(globals);

	std::cout << "Lagrangian Object Initialisation" << endl;
	lagrangian_object_vec[0].initialise(globals.PI);
	if (lagrangian_object_vec[0].type == 2) {
		tp.spring_network_output(globals.output_file, false, lagrangian_object_vec[0]);
	}else {
		tp.tecplot_output_lagrangian_object(lagrangian_object_vec[0], globals, domain, 0.0);
	}


    Solution soln;
    Solution residual;
	std::cout << "Importing Mesh" << endl;
     if(globals.mesh_type > 2){
        unstructured_mesh uns_mesh(domain,globals);
		std::cout << "Set Boundary Conditions" << endl;
        Boundary_Conditions bc(uns_mesh.get_num_bc());
        // assign external force terms
		std::cout << "Set External Forces" << endl;
        external_forces source_term(uns_mesh.get_total_cells());
        source_term.set_uniform_force(initial_conds.pressure_gradient); //pressure force applied through external force

		std::cout << "Set Initial Solution" << endl;
        //create solution
        soln.assign_memory(uns_mesh.get_total_cells());
        residual.assign_memory(uns_mesh.get_total_cells());
		std::cout << "Set Post Processor" << endl;
        post_processing post_processor(uns_mesh.get_total_cells());

		std::cout << "Initialising Solution" << endl;
        if(globals.restart_analysis == "true"){
            soln.import(globals);
			//lagrangian_object_vec[0].solution_read(globals);

        }else{
                soln.assign_pressure_gradient(initial_conds.rho_gradient,initial_conds.origin_loc,
                                    initial_conds.rho_origin_mag,uns_mesh,globals);
                soln.assign_velocity_gradient(initial_conds.vel_gradient,initial_conds.origin_loc,
                                        initial_conds.vel_origin_mag,uns_mesh,globals);
                soln.set_average_rho(initial_conds.average_rho);

        }

        bc.assign_boundary_conditions(uns_mesh,u_bcs);
		std::cout << "Running Solver" << endl;
		if (globals.gpu_solver == 1) {
			gpu.General_Purpose_Solver_mk_i(uns_mesh, soln, bc, source_term, globals, domain, initial_conds, u_bcs,
				mg, residual, fmg, post_processor, lagrangian_object_vec);

		}else {
			solve.General_Purpose_Solver_mk_i(uns_mesh,soln,bc,source_term,globals,domain,initial_conds,u_bcs,
										  mg,residual,fmg,post_processor);

		}

		//output object solution for restart
		if (lagrangian_object_vec[0].type == 2) {
			tp.spring_network_output(globals.output_file, false, lagrangian_object_vec[0]);
		}
		else {
			tp.tecplot_output_lagrangian_object(lagrangian_object_vec[0], globals, domain, 0.0);
		}

        soln.output(globals.output_file, globals, domain);
    //soln.output_centrelines(globals.output_file,globals,mesh);
    //post_processor.output_vtk_mesh(globals.output_file,globals,domain);

        std::clock_t end = clock();

        duration = double(end-start)/ CLOCKS_PER_SEC;

        output_globals(globals,duration,post_processor);

        std::cout << duration << std::endl;
        std::cout << "CPU Cycles" << double(end-start) << std::endl;

     }





}

void program::remove_existing_files(global_variables &globals){

    std::string output_location , temp;

    output_location = globals.output_file;


    temp = output_location + "/uy";
    remove_directory_files(temp);

    temp = output_location + "/vx";
     remove_directory_files(temp);

    temp = output_location + "/plt";
    remove_directory_files(temp);

    temp = output_location + "/plt/grid";
    remove_directory_files(temp);

	temp = output_location + "/plt/object";
	remove_directory_files(temp);

    }

void program::remove_directory_files(std::string output_location){

    fs::path p(output_location);

    std::string plt,dat,szplt;
    plt = ".plt";
    dat = ".dat";
    szplt = ".szplt";


    if(fs::exists(p) && fs::is_directory(p)){
        fs::directory_iterator end;

            for (fs::directory_iterator it(p); it !=end; ++it){
                try{
                    if( fs::is_regular_file(it->status()) && (it->path().extension().compare(plt) == 0 ||
                    it->path().extension().compare(dat) == 0 ||
                    it->path().extension().compare(szplt) == 0))
                    {
                        fs:remove(it->path());
                    }
                }catch(const std::exception &ex){
                    ex;
                }
            }
        }

    }

void program::output_globals (global_variables globals,double duration,post_processing &pp){

    std::string output_location;
    std::string filename;
    std::ofstream globals_txt ;
    std::string globals_file;
    output_location = globals.output_file;
    filename = globals.simulation_name;
    globals_file = output_location + "/globals.txt";
    double cycles;

    cycles = duration * CLOCKS_PER_SEC;
/// Generic Load Case Input

    globals_txt.open(globals_file.c_str(), ios::out);

    globals_txt << "Runtime:" << duration << "s" << endl;
    globals_txt << "CPU Cycles:" << cycles << endl;
    globals_txt << "File:" << filename.c_str()  << endl;

    globals_txt << "Tau:"  << globals.tau << endl;
    globals_txt << "Mach:" << globals.max_velocity *sqrt(3)<< endl;
    globals_txt << "Reynolds: " << globals.reynolds_number << endl;
    globals_txt << "Drag Coefficient: " << pp.drag_coef << endl;
    globals_txt << "Lift Coefficient: " << pp.lift_coef <<  endl;


    std::string reynolds_text;
    reynolds_text = globals.reynolds_number;

    globals_txt << "Separation Angles:" << endl;
    for(int i =0; i < pp.separation_angle.size(); i++){

        globals_txt << "grad magnitude: " << pp.separation_grad[i] << " separation angle: " <<
                pp.separation_angle[i] << " ," << pp.separation_i[i] << " ," << pp.separation_j[i]
                << "," << pp.separation_x[i] << "," << pp.separation_y[i] << endl;

    }

    globals_txt.close();


}


    // copy in binary mode
void program::copyfile( char* SRC,  std::string  DEST)
{
    DEST.append("/input.xml");

    std::ifstream src(SRC, std::ios::binary);
    std::ofstream dest(DEST, std::ios::binary);
    dest << src.rdbuf();

}


