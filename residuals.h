#ifndef RESIDUALS_H
#define RESIDUALS_H
#include "Solution.h"
#include "global_variables.h"
#include "cuda.h"
#include "cuda_runtime.h"

class residuals
{
    public:
        residuals();
        virtual ~residuals();
        void reset();
        void add_l2_norm_residuals(double f1, double rho, double f2, double u, double f3, double v);
        void add_l2_norm_residuals(Solution &residual,int i);
        void l2_norm_rms_moukallad(global_variables &globals);
		void l2_norm_rms_moukallad(global_variables &globals, double4 total_res);
		void l2_norm_rms_moukallad(global_variables &globals, double* rho, double* u, double* v, double* w, int n_blocks, int n_cells);
        void l2_norm_rms();
        double max_error();
		void add_l2_norm_residuals(double rho, double u, double v, double w, int i);

        void add_ansys_l2_norm_residuals(double f1, double rho, double f2, double u, double f3, double v,
                                          double delta_t);
        void ansys_5_iter_rms (int timestep);

		void force_totals(double* x, double* y, double* z, int n_blocks, int n_nodes);

		double force_x, force_y, force_z;  // object node forces
        double rho_rms, u_rms, v_rms, w_rms;
    protected:
    private:
     double rho_rms_num, u_rms_num, v_rms_num, w_rms_num, error;
    double rho_rms_den, u_rms_den, v_rms_den, w_rms_den;

    double max_rho, max_u, max_v;
};

#endif // RESIDUALS_H
