#include "Solution.h"
#include "gradients.h"
#ifndef POST_PROCESSING_H
#define POST_PROCESSING_H


class post_processing
{
    public:
        post_processing();
        virtual ~post_processing();
        post_processing(int _total_nodes);
        void calc_vorticity(Solution &x_grad, Solution &y_grad);
        void calc_streamfunction(Mesh &mesh,global_variables &globals,
            Boundary_Conditions &bcs);
            void Initialise();
        double *vorticity, *streamfunction;
        void cylinder_post_processing( unstructured_mesh &mesh, global_variables &globals,
                gradients &grads, Boundary_Conditions &bcs, Solution &soln,domain_geometry &geom, Solution &wall_shear_stress);
		void cylinder_post_processing(unstructured_mesh &mesh, global_variables &globals,
			gradients &grads, Boundary_Conditions &bcs, Solution &soln, domain_geometry &geom, double4 * res_face);
        double drag_coef, lift_coef;
        std::vector <double> separation_angle, separation_grad , separation_i, separation_j,separation_x,separation_y;

    protected:

    private:
         private:

        int total_nodes;
};

#endif // POST_PROCESSING_H
