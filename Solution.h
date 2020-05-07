#include "Mesh.h"
#include "vector_var.h"
#include <string>
#include "global_variables.h"
#include "domain_geometry.h"
#include "Boundary_Conditions.h"
#include "initial_conditions.h"
#include "unstructured_mesh.h"
#include "cuda.h"
#include "cuda_runtime.h"
#ifndef SOLUTION_H
#define SOLUTION_H


class Solution
{
    public:
        Solution(int _total_nodes);
        Solution();
        virtual ~Solution();
        void assign_memory(int _total_nodes);
        double get_rho( int i) {return rho[i];};
        double get_u( int i) {return u[i];};
       double get_v( int i) {return v[i];};
        double get_w( int i) {return w[i];};
        double get_average_rho (){return average_rho;};
       // double get_u_exact( int i) {return u_exact[i];};
       // double get_u_error( int i) {return error[i];};
        void set_average_rho(double arg) { average_rho = arg;};
        void set_rho( int i,double arg) {rho[i] =arg;};
        void set_u( int i,double arg) {u[i] =arg;};
        void set_v( int i,double arg) {v[i] =arg;};
        void set_w( int i,double arg) {w[i] =arg;};
        void assign_pressure_gradient( vector_var _gradient, vector_var gradient_origin,
            vector_var origin_magnitude,Mesh &Mesh, global_variables &globals);
        void assign_velocity_gradient( vector_var _gradient, vector_var gradient_origin,
            vector_var origin_magnitude,Mesh &Mesh,global_variables &globals);
        void uns_assign_pressure_gradient( vector_var _gradient, vector_var gradient_origin,
            vector_var origin_magnitude,unstructured_mesh &Mesh, global_variables &globals);
        void uns_assign_velocity_gradient( vector_var _gradient, vector_var gradient_origin,
            vector_var origin_magnitude,unstructured_mesh &Mesh,global_variables &globals);
        void update ( double rho, double u, double v, double w , int i);

        void output (std::string output_location, global_variables &globals,
        domain_geometry &geometry) ;
        void clone( Solution &soln_a);
		void clone(double4 *soln_a);
        void post_process(double pre_condition_gamma, Mesh &mesh, global_variables &globals,
                          initial_conditions &initials);
        void add_rho(int i, double arg) { rho[i] = rho[i] + arg;};
        void add_u(int i, double arg) { u[i] = u[i] + arg;}
        void add_v (int i , double arg) {v[i] = v[i] + arg;}
        void add_w (int i, double arg) {w[i] = w[i] + arg;}
        void restriction(Solution &soln, Mesh &coarse_mesh,
                           Mesh &fine_mesh, Boundary_Conditions &bc);
        void prolongation(Solution &coarse_soln, Solution &temp_soln, Solution &soln,
                            Mesh &coarse_mesh, Mesh &fine_mesh,
                            Boundary_Conditions &bc, bool fmg);
        void update_bcs(Boundary_Conditions &bcs,Mesh &mesh,domain_geometry &domain);
        void update_unstructured_bcs(Boundary_Conditions &bcs,unstructured_mesh &mesh,domain_geometry &domain,
                                    int time);
        void Initialise();
        void remove_double_errors();

        void update_gradients(Boundary_Conditions &bcs,Mesh &mesh,domain_geometry &domain,
                    int direction, Solution &src );
        void update_gradients_least_squares(Boundary_Conditions &bcs,Mesh &mesh,domain_geometry &domain,
                    int direction, Solution &src );
        void update_gradients_green_gauss(Boundary_Conditions &bcs,Mesh &mesh,domain_geometry &domain,
                    int direction, Solution &src );
        void output_centrelines (std::string output_location, global_variables &globals,
        Mesh &mesh,double time);
        void import(global_variables &globals);

		double *rho, *u, *v, *w;
		/*double *error, *u_exact;*/
		int total_nodes;
		double average_rho;


    protected:
    private:
       


};

#endif // SOLUTION_H
