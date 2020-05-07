#ifndef GLOBAL_VARIABLES_H
#define GLOBAL_VARIABLES_H
#include <string>
#include "domain_geometry.h"
#include <math.h>
#include "initial_conditions.h"

class global_variables
{
    public:
        global_variables();
        virtual ~global_variables();
        void initialise(domain_geometry domain, initial_conditions inits);
        void update_tau(domain_geometry domain);
        void update_visc(double Y );
        void update_fine_tau();
        void update_coarse_tau();
        void magnify_time_step();
        void reduce_time_step();
        std::string create_output_directory();

        //BC constants
        const int periodic  =3;
        const int dirichlet = 1;
        const int neumann = 2;

        //small number constant
        const double small_number = 0.0000000000000001;

        //PI
        const double PI = acos(-1.0);

        //tolerance in final solution
        double tolerance;

        double pre_conditioned_gamma = 1;
        double simulation_length = 2000; // cut off time in seconds
        double time_marching_step;
        double reynolds_number;
        double max_velocity;
        double knudsen_number;
        double tau;
        double arti_disp_kappa_2; // artiifical dissipation constant 2nd order
        double arti_disp_kappa_4;
        double martinelli;
        double ref_rho;
        double scale;
        double womersley_no;
        double visc;
	
        int max_mg_levels ; // number of coarse grids to enter
        int fmg_levels;
        int testcase; // couette = 1 , poiseuille = 2  // 5 is RBC inflation  /// 6 is Optical tweezers in fluid  /// 7 is optical tweezers (structural only)
						/// 8 is inviscid flow problems
        int mesh_type; // standard =1 , cosine = 2, unstructured_import =3 , unstructured_uniform =4 , cosine (GPU ) =5 , uniform (GPU) = 6
        int output_step;
        double fmg_tolerance;
		double fneq_min =0.0;
		int gpu_time_stepping; //local = 1, constant = 2, min =3, talts =4
		int gpu_solver; //1 = gpu, 0 = cpu
		int gpu_device = 0;

        std::string output_file_dir;
        std::string simulation_name = "";
        std::string output_file;
        std::string import_file;
        std::string import_format;
        std::string gradient_calc_type;
        const std::string test_case = "Couette";
        std::string  restart_analysis;
        std::string  time_stepping;
        std::string  preconditioning;
		bool preconditioning_active;

    protected:
    private:
};

#endif // GLOBAL_VARIABLES_H
