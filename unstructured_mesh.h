#ifndef UNSTRUCTURED_MESH_H
#define UNSTRUCTURED_MESH_H

#include "Mesh.h"
#include <vector>
#include "vector_var.h"
#include <string.h>
#include <tuple>

class unstructured_mesh : public Mesh
{
    public:
        unstructured_mesh(domain_geometry &domain ,global_variables &globals);
        virtual ~unstructured_mesh();
        int get_num_bc () { return num_bc_cells; };
        int get_n_vertices() { return n_vertices; };
        int get_n_cells () { return n_cells; };
         int get_n_faces () { return n_faces; };
         int get_n_neighbours () { return n_neighbours; };
		 int get_n_wall_cells() { return n_wall_cells; };
        int get_centre_node () { return centre_node; };
        double get_num_flux_calcs( int i) {return calcs_per_cell[i];};
        double get_flux_calc_directions( int i,int j) {return calculated_faces[i][j];};
         int get_bc_neighbour( int i) {return bc_neighbour[i];};
        double get_bc_face_x( int i) {return bc_face_x[i];};
		double get_delta_h() { return delta_h; };
        double get_bc_face_y( int i) {return bc_face_y[i];};
        double get_bc_face_z( int i) {return bc_face_z[i];};
        double get_bc_face_direction( int i) {return bc_face_direction[i];};
         std::vector <std::string> bc_types;
          int get_mesh_neighbour( int i) {return neighbour[i];};
            int get_mesh_owner( int i) {return owner[i];};
         double get_cell_cfl_x( int i) {return cell_cfl_x[i];};
        double get_cell_cfl_y( int i) {return cell_cfl_y[i];};
        double get_cell_cfl_z( int i) {return cell_cfl_z[i];}; 
        double get_cell_cfl_r( int i) {return cell_cfl_r[i];};
        std::vector <std::vector <int> > gradient_cells;
        std::vector <std::vector <int> > gradient_faces;
        double get_face_i( int i) {return face_i[i];};
        double get_face_j( int i) {return face_j[i];};
        double get_face_k( int i) {return face_k[i];};
        double get_face_x( int i) {return face_x[i];};
  double get_face_y( int i) {return face_y[i];};
        double get_face_z( int i) {return face_z[i];};
        double get_face_area( int i) {return face_area[i];};
        double get_delta_t_face( int i) {return delta_t_face[i];};

		double get_plateau_x(int i) { return plateau_x[i]; };
		double get_plateau_y(int i) { return plateau_y[i]; };
		double get_plateau_z(int i) { return plateau_z[i]; };

		int get_n_face_x() { return n_face_x; };
		int get_n_face_y() { return n_face_y; };
		int get_n_face_z() { return n_face_z; };
    protected:

    private:

         float * x, * y, * z;

		 double delta_h;
		 double * plateau_x, *plateau_y, *plateau_z;
         int * bc_neighbour;
         int * owner, * neighbour;
         std::vector <int> openfoam_map;
         std::vector <int> openfoam_map_B;

         double * bc_face_x,* bc_face_y, * bc_face_z, *bc_face_direction;

         double * face_x, *face_y,*face_z, *face_i,*face_j,*face_k,*face_area,*delta_t_face;
         double * cell_cfl_x, *cell_cfl_y, *cell_cfl_z, *cell_cfl_r;
         double * total_cell_area;
         int n_cells, n_vertices,  index_file,n_faces;
		 int n_face_x, n_face_y, n_face_z;
         int n_neighbours,n_face_points;
		 int n_wall_cells;
         int centre_node;
//        cgsize_t **ielem;
//         cgsize_t **ielem_mix;
//          cgsize_t **iparentdata;
//            cgsize_t ** i_elem_bc;

        std::vector<long long int> face_labels;
        std::vector <long long int>  ghost_faces;

        std::vector <std::vector <int> > calculated_faces ;
        std::vector <std::vector <int> > face_connectivity ;
        std::vector <std::vector <int> > volume_connectivity ;
		std::vector< std::tuple< double, double, int> > IBM_lookup;

        int * calcs_per_cell;


        int num_bc_cells;
        void calc_centre_node();
        void import_tecplot_mesh(global_variables &globals);
        void initialise_mesh_variables();
        void generate_internal_mesh(global_variables &globals);
        void destruct_mesh_variables(global_variables &globals);
        void import_cgns_mesh();
        void import_cgns_bcs();
        void calc_cell_center(std::vector<vector_var> &nodes, int i);
        void tet_volume (int a, int b, int c, int d, std::vector<vector_var> &nodes, int i);
        void tet_volume (int a, int b, int c, std::vector<vector_var> &nodes, int i);
        void calc_cell_volume(std::vector<vector_var> &nodes, int i);
        void tecplot_output_unstructured_grid(global_variables &globals);
        void tecplot_output_polyhedral_unstructured_grid(global_variables &globals);
        void calc_centroids_volume_area_normal(std::vector<vector_var> &nodes, int i);
        void cell_tet_area(std::vector<vector_var> &nodes, int i, int q, int r, int s,
        double &face_area, double &total_area,vector_var &centroid,vector_var &face_centroid
            ,vector_var &face_normal);
        void calc_face_centroid(std::vector<vector_var> &nodes, int i);
        void calc_delta_t(std::vector<vector_var> &nodes, int i,global_variables &globals);
        void calc_face_neighbours(global_variables &globals);
        void populate_face_neighbour(int index1, int index2);
        void get_num_bc_elements(global_variables &globals);
        void populate_ghost_cells(global_variables &globals);
        void open_cgns_file(global_variables &globals);
        void calc_unstructured_lattice_size();
        void face_neighbours_hash(int a, int b, int c, int d,int i,global_variables &globals);
        void get_ghost_face_hash(global_variables &globals);
        void populate_ghost_face_neighbour(int ghost, int index2);
         void read_bc_tags(global_variables &globals);
         void calculate_face_flux_calcs();
         void remove_face_from_flux_calcs(int i, int n);
         void output_mesh_to_text(global_variables &globals);
         void import_openfoam_mesh(global_variables &globals);
         void openfoam_connectivity(global_variables &globals);
         void import_openfoam_BC_tags(global_variables &globals);
         void populate_ghost_cell_from_foam_face(int face);
         void calc_internal_sphere_radius();
		 void create_uniform_openfoam_mesh(global_variables &globals,domain_geometry &domain);

		 void generate_immersed_boundary_look_up(global_variables &globals, domain_geometry &domain);


        void calc_face_delta_t(std::vector<vector_var> &nodes, int i,global_variables &globals);


		bool sort_tuple(const std::tuple <double, double, int>& a,
						const std::tuple <double, double, int>& b)
		{
			return (std::get<1>(a) > std::get<1>(b) && std::get<2>(a) > std::get<2>(b));
		}



         struct doCompare
           {
               doCompare( const unstructured_mesh& info ) : m_info(info) { } // only if you really need the object state
               const unstructured_mesh& m_info;

               bool operator()( const int & i, const int & j )
               {
                    double abs = 0.25;
                    // comparison code using m_info
                    if(fabs(m_info.x[i] -m_info.x[j]) < abs) {}
                    else if( m_info.x[i] < m_info.x[j] )return true;
                    else if((m_info.x[j] < m_info.x[i]) )return false;
                    else {};

                    if(fabs(m_info.y[i] -m_info.y[j]) < abs) {}
                    else if( m_info.y[i] < m_info.y[j] )return true;
                    else if((m_info.y[j] < m_info.y[i])  )return false;
                    else {};

                    if(fabs(m_info.z[i] -m_info.z[j]) < abs) {}
                    else if( m_info.z[i] < m_info.z[j] )return true;
                    else if((m_info.z[j] < m_info.z[i]) )return false;
                    else {};



                    if( (pow(m_info.x[i],2.0) + pow(m_info.y[i],2.0) + pow(m_info.z[i],2.0) ) <
                        (pow(m_info.x[j],2.0) + pow(m_info.y[j],2.0) + pow(m_info.z[j],2.0) )  ) return true;


                    return false;


               }
           };




};

#endif // UNSTRUCTURED_MESH_H
