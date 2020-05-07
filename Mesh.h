#ifndef Mesh_H
#define Mesh_H

#include "domain_geometry.h"
#include "global_variables.h"
class Mesh
{
    public:
        Mesh(domain_geometry domain, global_variables &globals);
        virtual ~Mesh();
        double get_node_x(int node_num);
        int get_num_x() { return num_x_cells; };
        int get_num_y () { return num_y_cells; };
        int get_total_cells () { return total_cells; };
        double get_dx() { return dx;};
        double get_dy() { return dy;};
        double get_X() { return X;};
        double get_Y() { return Y;};
        double get_Z() { return Z;};
        double get_centroid_x( int i) {return centroid_x[i];};
        double get_centroid_y( int i) {return centroid_y[i];};
        double get_centroid_z( int i) {return centroid_z[i];};
        double get_north_x( int i) {return north_x[i];};
        double get_north_y( int i) {return north_y[i];};
        double get_north_z( int i) {return north_z[i];};
        double get_east_x( int i) {return east_x[i];};
        double get_east_y( int i) {return east_y[i];};
        double get_east_z( int i) {return east_z[i];};
        double get_west_x( int i) {return west_x[i];};
        double get_west_y( int i) {return west_y[i];};
        double get_west_z( int i) {return west_z[i];};
        double get_south_x( int i) {return south_x[i];};
        double get_south_y( int i) {return south_y[i];};
        double get_south_z( int i) {return south_z[i];};
        double get_back_x( int i) {return back_x[i];};
        double get_back_y( int i) {return back_y[i];};
        double get_back_z( int i) {return back_z[i];};
        double get_front_x( int i) {return front_x[i];};
        double get_front_y( int i) {return front_y[i];};
        double get_front_z( int i) {return front_z[i];};
        double get_n_area( int i) {return n_area[i];};
        double get_s_area( int i) {return s_area[i];};
        double get_w_area( int i) {return w_area[i];};
        double get_e_area( int i) {return e_area[i];};
        double get_f_area( int i) {return f_area[i];};
        double get_b_area( int i) {return b_area[i];};
        double get_n_i( int i) {return n_i[i];};
        double get_n_j( int i) {return n_j[i];};
        double get_n_k( int i) {return n_k[i];};
        double get_s_i( int i) {return s_i[i];};
        double get_s_j( int i) {return s_j[i];};
        double get_s_k( int i) {return s_k[i];};
        double get_w_i( int i) {return w_i[i];};
        double get_w_j( int i) {return w_j[i];};
        double get_w_k( int i) {return w_k[i];};
        double get_e_i( int i) {return e_i[i];};
        double get_e_j( int i) {return e_j[i];};
        double get_e_k( int i) {return e_k[i];};
        double get_f_i( int i) {return f_i[i];};
        double get_f_j( int i) {return f_j[i];};
        double get_f_k( int i) {return f_k[i];};
        double get_b_i( int i) {return b_i[i];};
        double get_b_j( int i) {return b_j[i];};
        double get_b_k( int i) {return b_k[i];};
        double get_delta_t( int i) {return delta_t[i];};
        double get_delta_t_n( int i) {return delta_t_n[i];};
        double get_delta_t_e( int i) {return delta_t_e[i];};
        int get_n_node( int i) {return n_node[i];};
        int get_s_node( int i) {return s_node[i];};
        int get_e_node( int i) {return e_node[i];};
        int get_w_node( int i) {return w_node[i];};
        int get_f_node( int i) {return f_node[i];};
        int get_b_node( int i) {return b_node[i];};
        double get_cell_volume( int i) {return cell_volume[i];};

        void get_centroid(int i, vector_var &cell) ;

        domain_geometry create_coarse_mesh_domain( );


        void set_centroid_x( int i,double x) {centroid_x[i] = x;};
        void set_centroid_y( int i,double y) { centroid_y[i] = y;};
        void set_centroid_z( int i,double z) { centroid_z[i] = z;};

		double * cell_volume;

    protected:

         double * centroid_x, * centroid_y, * centroid_z;
        double * north_x, * north_y, * north_z;
        double * east_x, * east_y, * east_z;
        double * west_x, * west_y, * west_z;
        double * south_x, * south_y, * south_z;

        double * back_x, * back_y, * back_z;
        double * front_x, * front_y, * front_z;


        // magnitude of surface at cell interface
        double * n_area, * s_area, * w_area, * e_area;
        double * f_area, * b_area;

        

        double cosine_dt;
        //cell interface unit vectors
        double *n_i, *n_j, *n_k;
        double *s_i, *s_j, * s_k;
        double *w_i, *w_j, * w_k;
        double * e_i, *e_j, *e_k;
        double * f_i, *f_j, *f_k;
        double * b_i, *b_j, *b_k;
        // delta t max for each cell
        double * delta_t , *delta_t_n, *delta_t_e;

        //neighbouring nodes
        int * n_node, * s_node, * e_node, * w_node,* f_node, * b_node;
         void calc_unstructured_lattice_size();
        int total_cells;
          double X,Y,Z,dx,dy,dz,multi_grid_dt;
    private:


        void create_standard_mesh(int strt, int fnsh);
        void create_cosine_mesh(double PI,double dt);


        int num_x_cells, num_y_cells ;

        //centroid and cell interface locations


        // cell volume


        double cs;   //speed of sound in medium

};


#endif // Mesh_H
