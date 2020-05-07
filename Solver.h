#include "Mesh.h"
#include "Boundary_Conditions.h"
#include "Solution.h"
#include "vector_var.h"
#include "external_forces.h"
#include "global_variables.h"
#include "initial_conditions.h"
#include "flux_var.h"
#include "vector"
#include "post_processing.h"
#include "unstructured_bcs.h"
#include "gradients.h"
#ifndef SOLVER_H
#define SOLVER_H

using namespace std;

class Solver
{
    public:
        Solver();
        virtual ~Solver();
        void Uniform_Mesh_Solver( Mesh &Mesh , Solution &soln, Boundary_Conditions &bc,
                                   external_forces &source,global_variables &globals, domain_geometry &domain,
                                   initial_conditions &init_conds,quad_bcs_plus &quad_bcs_orig, int mg,
                                   Solution &residual,int fmg, post_processing &pp);
        void Unstructured_Mesh_Solver( Mesh &Mesh , Solution &soln, Boundary_Conditions &bc,
                                   external_forces &source,global_variables &globals, domain_geometry &domain,
                                   initial_conditions &init_conds,quad_bcs_plus &quad_bcs_orig, int mg,
                                   Solution &residual,int fmg, post_processing &pp);
        void Unstructured_Mesh_Solver_mk_ii( Mesh &Mesh , Solution &soln, Boundary_Conditions &bc,
               external_forces &source,global_variables &globals, domain_geometry &domain,
               initial_conditions &init_conds,quad_bcs_plus &quad_bcs_orig, int mg,
               Solution &residual,int fmg, post_processing &pp);
         void General_Purpose_Solver_mk_i( unstructured_mesh &Mesh , Solution &soln, Boundary_Conditions &bc,
               external_forces &source,global_variables &globals, domain_geometry &domain,
               initial_conditions &init_conds,unstructured_bcs &quad_bcs_orig, int mg,
               Solution &residual,int fmg, post_processing &pp);


        void multi_grid_agglomoration( Solution &residuals , Solution &soln,
                                         int cycle_no, Mesh &fine_mesh,  quad_bcs_plus &bcs,
                                         initial_conditions &init_conds, int &mg, global_variables globals,
                                         domain_geometry &fine_domain,Boundary_Conditions &fine_bc
                                         , post_processing &pp);
        void populate_e_alpha(std::vector<vector_var> &e_alpha,double * lattice_weight, double c, double PI, int k );
        void inverse_weighted_distance_interpolation(double &u, double &v, double &rho, Boundary_Conditions &bcs,
                                Mesh &Mesh , domain_geometry &domain,Solution &soln
                                ,vector_var &cell_interface, int k,int i, int neighbour,
                                std::vector<vector_var> &e_alpha, int j,std::vector<int> &cell_nodes);
        void get_cell_nodes(std::vector<int> &cell_nodes, Boundary_Conditions &bcs,int neighbour,
                                Mesh &Mesh,int i, int j );
        void get_cfl(double &delta_t,Solution &soln,unstructured_mesh &Mesh, global_variables &globals,
                        double* delta_t_local, int* delta_t_frequency, Solution &cfl_areas);

        void find_real_time( double* delta_t_local, double* local_time, bool* calc_face,
                    unstructured_mesh &Mesh,bool* calc_cell);

    // initialise variables
    protected:
    private:

        double dt;
        double tau;
        double kine_viscosity;
        double c,cs;

//        struct vector_var {
//            double x;
//            double y;
//            double z;
//        };
        vector_var get_e_alpha(int k, double &lattice_weight, double c,double PI );

        void get_cell_gradients(Mesh &Mesh, int i, int neighbour, int j, Solution &temp_soln,
                                vector_var &delta_rho, vector_var &delta_rho1,
                                vector_var &delta_u, vector_var &delta1_u1,
                                vector_var &delta_v, vector_var &delta_v1,
                                Boundary_Conditions &bcs);


        struct bc_var{

            bool present;
            double rho;
            double u;
            double v;
            int vel_type;
            int rho_type;
            int periodic_node;

        };



        void truncate_flux(flux_var &flux);
        void truncate_flux(double &val);

        void get_gradients_weighted_average( gradients &grads, int i, int neighbour, double m1,
                double m2, vector_var &u,vector_var &v,vector_var &w,vector_var &rho);

         void cell_interface_initialiser( double &rho_interface,vector_var &rho_u_interface,
                                        flux_var &x_flux,flux_var &y_flux);

        double feq_calc(double weight, vector_var e_alpha, vector_var u_lattice, double u_magnitude,
                        double cs, double rho_lattice);

        double feq_calc_incomp(double weight, vector_var e_alpha, vector_var u_lattice, double u_magnitude,
                        double cs, double rho_lattice, double rho_0,int k);

        void cell_interface_variables( int j, int i, vector_var &interface_node, int &neighbour, double &interface_area,
                              vector_var &cell_normal, Boundary_Conditions &boundary_conditions,  bc_var &bc,
                              Mesh &Mesh, vector_var &cell_2);
        void populate_feq(double u_lattice[],double v_lattice[],
            double w_lattice[],double rho_lattice[],double lattice_weight[],
            double feq_lattice[],int k,global_variables &globals);

        void get_weighted_average( gradients &grads, int i, int neighbour, double m1, double m2,
            vector_var &u, vector_var &v, vector_var &w, vector_var &rho, unstructured_mesh &mesh);

        void populate_lattice_macros(double u_lattice[],double v_lattice[],
            double w_lattice[],double rho_lattice[],vector_var cell_1, vector_var cell_2,
            vector_var interface_node, int i, int neighbour, gradients &grads,Solution &temp_soln,
           unstructured_mesh &mesh, Boundary_Conditions &bcs,vector_var &cell_normal);

        void calculate_flux_at_interface(vector_var u_interface, double dt, global_variables &globals,
                double local_viscosity[], double rho_interface,double lattice_weight[],
                flux_var &cell_flux,int i,double feq_lattice[],vector_var cell_normal,double interface_area, double local_fneq[]);

        void cell_interface_variables( int face, int i, vector_var &interface_node, int &neighbour, double &interface_area,
                              vector_var &cell_normal, Boundary_Conditions &boundary_conditions,  bc_var &bc,
                              unstructured_mesh &Mesh, vector_var &cell_2, vector_var &cell_1 );
		void populate_cfl_areas(Solution &cfl_areas, unstructured_mesh &Mesh);
};

#endif // SOLVER_H
