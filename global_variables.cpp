#include "global_variables.h"
#include "domain_geometry.h"
#include <boost/algorithm/string/replace.hpp>
#include <boost/filesystem.hpp>
#include <sstream>

using namespace boost::filesystem;
using namespace std;

global_variables::global_variables()
{
    //ctor
}

global_variables::~global_variables()
{
    //dtor
}

void global_variables::initialise(domain_geometry domain,initial_conditions initial_conds){
    //tau = 0.5 + viscosity / dt/ gamma
    // viscosity = MA/root(3) /Re
     std::ostringstream s;
     visc =  max_velocity * domain.ref_length/reynolds_number;

     tau = 3*visc/domain.dt + 0.5;  // non-dimensional dt is 1 here
//    tau = 0.5 + initial_conds.average_rho*max_velocity*3/reynolds_number *
//                domain.Y/domain.dt*pre_conditioned_gamma;
    knudsen_number = max_velocity *sqrt(3) / reynolds_number;

	if (testcase == 8) {
		tau = 0.5;
		visc = 0.0;
	}


    s << "RE_" << reynolds_number << " N_CELLS_Y" << domain.Y << " N_CELLS_x" << domain.X <<
                    " MA_" << max_velocity *sqrt(3)/scale << " dt_" << domain.dt
                    << " DT_" << time_marching_step << " wom_ " << womersley_no;
    simulation_name = s.str();
    boost::replace_all(simulation_name,".","_");
    output_file = create_output_directory();
    simulation_length = simulation_length *scale;
    ref_rho = initial_conds.average_rho;

    if( preconditioning == "true"){
		preconditioning_active = true;
		if (fabs(pre_conditioned_gamma - 1) < 0.00001) {
			pre_conditioned_gamma = max_velocity * 2.0 *sqrt(3);
		}

        
	}
	else {
		preconditioning_active = false;

	}

	if (time_stepping == "local") {
		gpu_time_stepping = 1;
	}
	else if (time_stepping == "constant") {
		gpu_time_stepping = 2;
	}
	else if (time_stepping == "min") {
		gpu_time_stepping = 3;
	}



}
void global_variables::update_visc(double y){

    visc =  max_velocity * y/reynolds_number;
}
void global_variables::update_coarse_tau(){

    tau = tau/2.0;
}

void global_variables::update_fine_tau(){

    tau = tau*2.0;
}

void global_variables::update_tau( domain_geometry domain){
     tau = 0.5 + max_velocity/reynolds_number * ceil(domain.Y /domain.dt) *pre_conditioned_gamma;

}
void global_variables::magnify_time_step( ){
     time_marching_step = time_marching_step * 2;

}

void global_variables::reduce_time_step( ){
     time_marching_step = time_marching_step / 2;

}
std::string global_variables::create_output_directory(){

    std::string output_file,plt,ux,vy,max_u, object;
    std::string folder;
    std::ostringstream s ,s1,s2,s3,s4,s5;
    //output_file = "C:/Users/brendan/Dropbox/PhD/Test Cases/Poiseuille Flow/";

    output_file = output_file_dir;

    if( simulation_name ==  "default"){
        s << "tol " << tolerance << " RE " << reynolds_number
         << " t " << time_marching_step << " gamma " << pre_conditioned_gamma << " tau "
         << tau << " wom " << womersley_no;
         folder = s.str();

    }else{

        folder = simulation_name;

    }

    //folder.replace(folder.begin(),folder.end(), ".",  "_");
    boost::system::error_code ec;
    boost::replace_all(folder, "." , "_");
    output_file = output_file + folder;

    boost::filesystem::path dir(output_file);
    boost::filesystem::create_directories(dir);

    // create plt, uy, vx folders
    s1 << "/plt";
    plt = output_file + s1.str();
    boost::filesystem::path dir1(plt);
    boost::filesystem::create_directories(dir1);

	s5<< "/plt/object";
	object = output_file + s5.str();
	boost::filesystem::path dir5(object);
	boost::filesystem::create_directories(dir5);

    s1 << "/grid";
    plt = output_file + s1.str();
    boost::filesystem::path dir4(plt);
    boost::filesystem::create_directories(dir4);

    // create plt, uy, vx folders
    s2 << "/uy";
    ux = output_file + s2.str();
    boost::filesystem::path dir2(ux);
    boost::filesystem::create_directories(dir2);

    // create plt, uy, vx folders
    s3 << "/vx";
    vy = output_file + s3.str();
    boost::filesystem::path dir3(vy);
    boost::filesystem::create_directories(dir3);

    output_file = output_file;
    return output_file;

}
