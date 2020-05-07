#include "Solution.h"
#include "Boundary_Conditions.h"
#include "flux_var.h"
#include "global_variables.h"
#ifndef ARTIFICIAL_DISSIPATION_H
#define ARTIFICIAL_DISSIPATION_H


class artificial_dissipation
{
    public:
        artificial_dissipation();
        artificial_dissipation(int nodes,global_variables &globals);
        virtual ~artificial_dissipation();
        double * global_JST_switch_x,* global_JST_switch_y;
        void get_global_jst(Solution &soln, Boundary_Conditions &bcs,
                            Mesh &Mesh, domain_geometry &domain);
        void set_local_jst();
        void reset_local_jst_switch();
        double local_jst_switch_x, local_jst_switch_y;
        double kappa_2, kappa_4;
        void get_local_coeffs(Solution &soln, Boundary_Conditions &bcs,
                                            Mesh &Mesh, Solution &local_soln,domain_geometry &domain,
                                            int j, int i);

        double global_2nd_order,global_4th_order;
        double spectral_radii[4] ; // 1 and 2 = i+ 1;  3 + 4 = j+1
        double martinelli_exponent;
        flux_var local_flux;

    protected:

    private:
        double jst_num,jst_den;

        double maximum( double m1, double zero, double p1, double p2);
        double second_order_difference (int var, int neighbour, int i,Solution &soln, Solution &local_soln);
        double _4th_order_difference(int var, int neighbour, int i,Solution &soln, Boundary_Conditions &bcs,
                                            Mesh &Mesh, Solution &local_soln, domain_geometry &domain,int j);

};

#endif // ARTIFICIAL_DISSIPATION_H
