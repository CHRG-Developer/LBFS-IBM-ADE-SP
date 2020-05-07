#include "residuals.h"
#include <math.h>

#include "global_variables.h"
residuals::residuals()
{
    //ctor
    max_rho = 0;
    max_u = 0;
    max_v = 0;
}

residuals::~residuals()
{
    //dtor
}

void residuals::reset(){

        rho_rms = 0;
        u_rms =0;
        v_rms = 0;
        w_rms = 0;
        rho_rms_num =0 ;
        u_rms_num =0;
        v_rms_num = 0;
        w_rms_num = 0;
        rho_rms_den =0;
        u_rms_den = 0;
        v_rms_den = 0;
        w_rms_den = 0;

}

void residuals::add_l2_norm_residuals(Solution &residual,int i){

    rho_rms_num = rho_rms_num + pow((residual.get_rho(i)),2.0);
    u_rms_num = u_rms_num + pow((residual.get_u(i)),2.0);
    v_rms_num = v_rms_num + pow((residual.get_v(i)),2.0);
    w_rms_num = w_rms_num + pow((residual.get_w(i)),2.0);
    rho_rms_den = rho_rms_den + 1; // get number of cells


}

void residuals::add_l2_norm_residuals(double rho, double u, double v, double w, int i) {

	rho_rms_num = rho_rms_num + rho* rho;
	u_rms_num = u_rms_num + u*u;
	v_rms_num = v_rms_num + v*v;
	w_rms_num = w_rms_num + w*w;
	rho_rms_den = rho_rms_den + 1; // get number of cells


}


void residuals::add_l2_norm_residuals(double f1, double rho, double f2, double u,
        double f3, double v){

            rho_rms_num = rho_rms_num + pow( f1 -rho ,2.0 );
            if( f1 < pow(1,-8)){
                rho_rms_den = rho_rms_den + pow(1,5) ;

            }else{
                rho_rms_den = rho_rms_den + pow(f1,2) ;

            }
            u_rms_num = u_rms_num + pow( f2 -u,2.0 );
            u_rms_den = u_rms_den + pow(f2,2) ;

            v_rms_num = v_rms_num + pow( f3 -v ,2.0 ) ;
            v_rms_den = v_rms_den +pow(f3,2);
            v_rms_den = v_rms_den + 1;

    };


void residuals::l2_norm_rms_moukallad(global_variables &globals){

        double ref_rho = globals.ref_rho;
        double ref_mom = globals.ref_rho * globals.max_velocity;
        double n_cells = rho_rms_den;

         rho_rms = sqrt(rho_rms_num / n_cells) /ref_rho ;

         u_rms = sqrt(u_rms_num/n_cells) /ref_mom;
         v_rms = sqrt(v_rms_num/n_cells) /ref_mom;
         w_rms = sqrt(w_rms_num/n_cells) /ref_mom;

}

void residuals::force_totals( double* x, double* y, double* z,  int n_blocks, int n_nodes) {

	force_x = 0.0;
	force_y = 0.0;
	force_z = 0.0;

	//sum blocksums
	for (int i = 0; i < n_blocks; i++) {
		force_x = force_x + x[i];
		force_z = force_z + z[i];
		force_y = force_y + y[i];
	}
	force_x = sqrt(force_x);
	force_y = sqrt(force_y);
	force_z = sqrt(force_z);
	
	return;
}



void residuals::l2_norm_rms_moukallad(global_variables &globals, double* rho, double* u, double* v, double* w, int n_blocks,int n_cells) {

	double ref_rho = globals.ref_rho;
	double ref_mom = globals.ref_rho * globals.max_velocity;
	
	//sum blocksums
	for (int i = 0; i < n_blocks; i++) {
		rho_rms = rho_rms + rho[i];
		u_rms = u_rms + u[i];
		v_rms = v_rms + v[i];
		w_rms = w_rms + w[i];
	}
	//get rms
	rho_rms = sqrt(rho_rms / n_cells) / ref_rho;
	u_rms = sqrt(u_rms / n_cells) / ref_mom;
	v_rms = sqrt(v_rms / n_cells) / ref_mom;
	w_rms = sqrt(w_rms / n_cells) / ref_mom;

}


void residuals::l2_norm_rms(){

        global_variables globals;
        if (rho_rms_den < globals.small_number){
            rho_rms = 0 ;
        }else{
            rho_rms = sqrt(rho_rms_num / rho_rms_den);
        }

        if (u_rms_den < globals.small_number){
            u_rms = u_rms ;
        }else{
            u_rms = sqrt(u_rms_num / u_rms_den);
        }

        if (v_rms_den < globals.small_number){
            v_rms = v_rms;
        }else{
            v_rms = sqrt(v_rms_num / v_rms_den) ;
        }

}

double residuals::max_error (){

    double temp;
	temp = w_rms;
    temp = fmax(w_rms, u_rms);
    
    return fmax(temp, v_rms);

}

void residuals::add_ansys_l2_norm_residuals(double f1, double rho, double f2, double u, double f3,
 double v,double factor){

    rho_rms_num = rho_rms_num + pow( (f1 -rho)/factor ,2.0 );
    u_rms_num = u_rms_num + pow( (f2 - u)/factor ,2.0 );
    v_rms_num = v_rms_num + pow( (f3 - v)/factor ,2.0 ) ;


};


// get max residual from first 5 timesteps and use this to scale subsequent timesteps
void residuals::ansys_5_iter_rms(int timestep){

    if (timestep < 5){

        max_rho = fmax(rho_rms_num,max_rho);
        max_u = fmax(u_rms_num, max_u);
        max_v = fmax(v_rms_num, max_v);

    }

     global_variables globals;

    if (max_rho < globals.small_number){
        rho_rms = rho_rms ;
    }else{
        rho_rms = sqrt(rho_rms_num / max_rho);
    }

    if (max_u < globals.small_number){
        u_rms = u_rms ;
    }else{
        u_rms = sqrt(u_rms_num / max_u);
    }

    if (max_v < globals.small_number){
        v_rms = v_rms;
    }else{
        v_rms = sqrt(v_rms_num / max_v) ;
    }



}
