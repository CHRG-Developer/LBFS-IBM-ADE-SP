#include "quad_bcs.h"
#include "quad_bcs_plus.h"
#include "unstructured_mesh.h"
#include "unstructured_bcs.h"

#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H


class Boundary_Conditions
{
    public:
        Boundary_Conditions(int total_cells);
        virtual ~Boundary_Conditions();
        void assign_boundary_conditions(int num_x, int num_y, quad_bcs);
        void assign_boundary_conditions(int num_x, int num_y, quad_bcs_plus _bc,
            int testcase);
       void assign_boundary_conditions(unstructured_mesh &mesh, unstructured_bcs _bc);

        bool get_bc( int i) {return bc[i];};
        bool get_bc_include( int i) {return bc_include[i];};


        double get_rho( int i) {return rho[i];};
      double get_u( int i) {return u[i];};
        double get_v ( int i) {return v [i];};
        double get_w ( int i) {return w[i];};
       

        int get_rho_type( int i) {return type_rho[i];};

        int get_vel_type( int i) {return type_vel[i];};

        int get_periodic_node( int i) {return periodic_node[i];};
       int get_neighbour(int i){return neighbour[i];};
       std::string get_name( int i) {return name[i];};
    protected:
    private:

        //centroid and cell interface locations

        bool *bc,*bc_include;
        double * rho, * u, * v , * w;
        int *type_rho;
        int *type_vel;
        int *periodic_node;
        int *neighbour;

       
        std::string *name;




        };

#endif // BOUNDARY_CONDITIONS_H
