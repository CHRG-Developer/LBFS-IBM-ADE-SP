#ifndef GRADIENTS_H
#define GRADIENTS_H

#include "Mesh.h"
#include "vector_var.h"
#include <string>
#include "global_variables.h"
#include "domain_geometry.h"
#include "Boundary_Conditions.h"
#include "initial_conditions.h"
#include "unstructured_mesh.h"
#include "Solution.h"




class gradients
{
    public:
        gradients();
        gradients(int total_cells);
        virtual ~gradients();
        void Initialise();
        void Get_LS_Gradients(Boundary_Conditions &bcs,unstructured_mesh &mesh,domain_geometry &domain,
                      Solution &src, global_variables &globals );
        void add_LS_contributions(vector_var &RHS_rho, vector_var &RHS_u, vector_var &RHS_v, vector_var &RHS_W, Boundary_Conditions &bcs,unstructured_mesh &mesh,domain_geometry &domain,
                      Solution &src , int i1, int i,int nb);
        void get_ghost_grads( unstructured_mesh &mesh,Solution &src , int ghost, int i,bool rho_type,
                    Boundary_Conditions &bcs,int face,global_variables &globals);
        void get_ghost_grads_limits( unstructured_mesh &mesh,Solution &src , int ghost, int i,bool rho_type,
                    Boundary_Conditions &bcs,int face,global_variables &globals);

        void add_non_macrovariable_contributions( Boundary_Conditions &bcs,
                        unstructured_mesh &mesh,domain_geometry &domain,
                      Solution &src , int i1, int i,double beta_0, int face,int i_nb) ;
        void pre_fill_LHS_and_RHS_matrix(Boundary_Conditions &bcs,unstructured_mesh &mesh,domain_geometry &domain,
                      Solution &src, global_variables &globals );


        double limit_grad(double min_u, double max_u, double grad_u );
        vector_var get_u(int i) { return u[i];};
        vector_var get_v(int i) { return v[i];};
        vector_var get_w(int i) { return w[i];};
        vector_var get_rho(int i) { return rho[i];};

		vector_var * u, *v, *w, *rho;
		double *LHS_xx, *LHS_xy, *LHS_xz, *LHS_yx, *LHS_yy, *LHS_yz, *LHS_zx, *LHS_zy, *LHS_zz;
		double *RHS_x, *RHS_y, *RHS_z;

    protected:

    private:
     





     int total_nodes;
};

#endif // GRADIENTS_H
