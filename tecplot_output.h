#ifndef TECPLOT_OUTPUT_H
#define TECPLOT_OUTPUT_H
#include "global_variables.h"
#include "Mesh.h"
#include "Solution.h"
#include "Boundary_Conditions.h"
#include "post_processing.h"
#include "lagrangian_object.h"

class tecplot_output
{
    public:
        tecplot_output();
        tecplot_output(global_variables &globals, Mesh &Mesh, Solution &Soln,
                             Boundary_Conditions &bcs, int fileType_e, double timestamp,
                             post_processing &pp);
		
        void tecplot_output_unstructured_soln(global_variables &globals, unstructured_mesh &Mesh, Solution &Soln,
                             Boundary_Conditions &bcs,  int time,
                             post_processing &pp, Solution &residual, double * local_delta_t, double * local_fneq);
        virtual ~tecplot_output();

		void tecplot_output_lagrangian_object(lagrangian_object &object, global_variables &globals, domain_geometry &geometry,double timestamp);
		void tecplot_output_lagrangian_object_gpu(double * vel_x, double * vel_y, double * vel_z, double * x, double * y, double * z, double * force_x, double * force_y, double * force_z,
			global_variables &globals, domain_geometry &geometry, double timestamp, std::string obj_name, int obj_nodes, int depth_nodes, int radial_nodes);


		
		void spring_network_output(std::string fileName, bool Ctn, lagrangian_object &object);
		
		void spring_network_output_gpu(std::string fileName, bool Ctn, double * obj_vel_x, double * obj_vel_y, double * obj_vel_z, double * obj_x, double * obj_y, double * obj_z,
			double * obj_force_x, double * obj_force_y, double * obj_force_z, std::string name, int time, int num_nodes, int num_tets, int * tet_connectivity, double * nodal_area,bool final, double * curvature);
    protected:

    private:
};

#endif // TECPLOT_OUTPUT_H
