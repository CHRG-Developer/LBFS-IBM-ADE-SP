#include "artificial_dissipation.h"
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>


artificial_dissipation::artificial_dissipation()
{
    //ctor
}
artificial_dissipation::artificial_dissipation(int total_nodes, global_variables &globals)
{
    //ctor
     global_JST_switch_x = new double [total_nodes +1];
        if (global_JST_switch_x==NULL) exit (1);
    global_JST_switch_y = new double [total_nodes +1];
    if (global_JST_switch_y==NULL) exit (1);

    martinelli_exponent = globals.martinelli;
    kappa_2 = globals.arti_disp_kappa_2;
    kappa_4 = globals.arti_disp_kappa_4;

}

artificial_dissipation::~artificial_dissipation()
{
    //dtor
    delete [] (global_JST_switch_x);
    global_JST_switch_x = NULL;
    delete [] (global_JST_switch_y);
    global_JST_switch_y = NULL;

}

void artificial_dissipation::get_global_jst(Solution &soln, Boundary_Conditions &bcs,
                                            Mesh &Mesh, domain_geometry &domain)
{
    int neighbour;

    for(int i = 0; i < Mesh.get_total_cells(); i++){

        if (bcs.get_bc(i)){
            global_JST_switch_x[i] = 0.0;
             global_JST_switch_y[i] = 0.0;

        }else{


        //first go x-direction
        // m1 = minus 1
        //p1 = plus1
        // zero = current node
        double m1,p1,zero;
        neighbour = Mesh.get_w_node(i);
        zero = soln.get_rho(i);
        m1 = soln.get_rho(neighbour);
        neighbour = Mesh.get_e_node(i);
        p1 = soln.get_rho(neighbour);

        global_JST_switch_x[i] = fabs( (m1 - 2*zero + p1)/ (m1 + 2*zero + p1));
        neighbour = Mesh.get_s_node(i);
        m1 = soln.get_rho(neighbour);

        neighbour = Mesh.get_n_node(i);
        p1 = soln.get_rho(neighbour);



        global_JST_switch_y[i] = fabs( (m1 - 2*zero + p1)/ (m1 + 2*zero + p1));
        }
    }



}
void artificial_dissipation::get_local_coeffs(Solution &soln, Boundary_Conditions &bcs,
                                            Mesh &Mesh, Solution &local_soln, domain_geometry &domain,
                                            int j, int i)
{
    int neighbour;
    double phi_i, phi_i_p1;
    double lambda_flux;

        if (bcs.get_bc(i)){
            global_2nd_order = 0.0;
             global_4th_order = 0.0;
             local_flux.P =0.0;
             local_flux.momentum_x =0.0;
             local_flux.momentum_y = 0.0;
             local_flux.momentum_z = 0.0;


        }else{

        /// Get dissipation coefficients
        double m1,p1,zero,p2;
        jst_num = -2*local_soln.get_rho(i);
        jst_den = 2*local_soln.get_rho(i);
        if( j == 0 || j == 2){


            neighbour = Mesh.get_w_node(i);
            m1 = global_JST_switch_x[neighbour];
            jst_num = jst_num + soln.get_rho(neighbour);
            jst_den = jst_den  + soln.get_rho(neighbour);
            neighbour = Mesh.get_e_node(i);
            p1 = global_JST_switch_x[neighbour];
            jst_num = jst_num + soln.get_rho(neighbour);
            jst_den = jst_den + soln.get_rho(neighbour);

            neighbour = Mesh.get_e_node(neighbour);
            p2 = global_JST_switch_x[neighbour];

        }else{

            neighbour = Mesh.get_s_node(i);
            m1 = global_JST_switch_y[neighbour];
            jst_num = jst_num + soln.get_rho(neighbour);
            jst_den = jst_den  + soln.get_rho(neighbour);
            neighbour = Mesh.get_n_node(i);
            p1 = global_JST_switch_y[neighbour];
            jst_num = jst_num + soln.get_rho(neighbour);
            jst_den = jst_den + soln.get_rho(neighbour);
            neighbour = Mesh.get_n_node(neighbour);
            p2 = global_JST_switch_y[neighbour];

        }

        zero = jst_num/jst_den;

        global_2nd_order = std::max(m1,std::max(zero,std::max(p1,p2))) * kappa_2;
        global_4th_order = std::max(0.0,(kappa_4 - global_2nd_order));




        /// Get spectral radii for euler equations (inviscid)

        /// needs allowances for pre conditioning


        /// artificial dissipation calcs to get lambda x

         // spectral radii for cell center
        spectral_radii[0] = fabs( local_soln.get_u(i))+ domain.cs;  // x-direction
        spectral_radii[2] = fabs( local_soln.get_v(i))+ domain.cs;  // y-direction




        switch(j){

            case 0: //west
                neighbour = Mesh.get_w_node(i);
                break;
            case 1: // south
                neighbour = Mesh.get_s_node(i);
                break;
            case 2: // east
                neighbour = Mesh.get_e_node(i);
                break;
            case 3: // north
                neighbour = Mesh.get_n_node(i);
                break;

            }

        //neighbour node spectral radii
        spectral_radii[1] = fabs( soln.get_u(neighbour))+ domain.cs; // x direction
        spectral_radii[3] = fabs( soln.get_v(neighbour))+ domain.cs; // y-direction



        // Martinelli scaling

        switch(j){

            case 0:
            case 2: //west or east
                    phi_i = 1 +pow( (spectral_radii[2]/spectral_radii[0]), martinelli_exponent) ;

                    phi_i_p1 = 1 + pow((spectral_radii[3]/spectral_radii[1]),martinelli_exponent);

                    lambda_flux = 0.5* ( phi_i * spectral_radii[0] + phi_i_p1 *spectral_radii[1]);

                break;
            case 1:
            case 3: // south or north
                    phi_i = 1 +pow( (spectral_radii[0]/spectral_radii[2]), martinelli_exponent) ;

                    phi_i_p1 = 1 + pow((spectral_radii[1]/spectral_radii[3]),martinelli_exponent);
                     lambda_flux = 0.5* ( phi_i * spectral_radii[2] + phi_i_p1 *spectral_radii[3]);
                break;


            }


        local_flux.P =  lambda_flux *global_2nd_order
                    * second_order_difference(1,neighbour,i,soln,local_soln)
                - lambda_flux  * global_4th_order
                * _4th_order_difference(1,neighbour,i,soln,bcs,Mesh,local_soln,domain,j);
        local_flux.momentum_x =  lambda_flux *global_2nd_order
                * second_order_difference(2,neighbour,i,soln,local_soln)
                - lambda_flux  * global_4th_order
                * _4th_order_difference(2,neighbour,i,soln,bcs,Mesh,local_soln,domain,j);
        local_flux.momentum_y = lambda_flux *global_2nd_order
                * second_order_difference(3,neighbour,i,soln,local_soln)
                - lambda_flux  * global_4th_order
                * _4th_order_difference(3,neighbour,i,soln,bcs,Mesh,local_soln,domain,j);
        local_flux.momentum_z = 0.0;




    }



}



void artificial_dissipation::reset_local_jst_switch(){
    local_jst_switch_x = 0.0;
    local_jst_switch_y = 0.0;
    jst_num = 0.0;
    jst_den = 0.0;

    local_flux.P = 0.0;
    local_flux.momentum_x = 0.0;
    local_flux.momentum_y = 0.0;
    local_flux.momentum_z = 0.0;
}


//use local_soln for current volume as this is the must up to date in the multi stage solution
double artificial_dissipation::second_order_difference(int var, int neighbour, int i,Solution &soln,
                                                       Solution &local_soln){

                    double temp;

                    switch(var){

                        case(1):
                            temp = soln.get_rho(neighbour) - local_soln.get_rho(i);
                            break;
                        case(2):
                            temp = soln.get_u(neighbour) - local_soln.get_u(i);
                            break;
                        case(3):
                            temp = soln.get_v(neighbour) - local_soln.get_v(i);
                            break;
                        case(4):
                            temp = soln.get_w(neighbour) - local_soln.get_w(i);
                            break;


                    }

                           return temp;

}

double artificial_dissipation::_4th_order_difference(int var, int neighbour, int i,Solution &soln,
                        Boundary_Conditions &bcs,Mesh &Mesh, Solution &local_soln,
                        domain_geometry &domain,int j){

                    double temp;

                    switch(var){

                        case(1):

                            switch(j){
                                case(0):
                                    if(bcs.get_bc(neighbour)){
                                            temp = 2*local_soln.get_rho(i)      // U j
                                                -soln.get_rho(neighbour)      // U j+1
                                                - soln.get_rho(Mesh.get_e_node(neighbour));

                                    }else{
                                        temp = soln.get_rho(Mesh.get_w_node(neighbour) ) // U j+2
                                        - 3* soln.get_rho(neighbour)                //U j+1
                                        + 3* local_soln.get_rho(i)                  // U j
                                        - soln.get_rho(Mesh.get_e_node(neighbour));  // U j-1

                                    }

                                    break;
                                case(1):
                                      if(bcs.get_bc(neighbour)){
                                            temp = 2*local_soln.get_rho(i)      // U j
                                                -soln.get_rho(neighbour)      // U j+1
                                                - soln.get_rho(Mesh.get_n_node(neighbour));

                                        }else{

                                        temp = soln.get_rho(Mesh.get_s_node(neighbour) ) // U j+2
                                            - 3* soln.get_rho(neighbour)                //U j+1
                                            + 3* local_soln.get_rho(i)                  // U j
                                            - soln.get_rho(Mesh.get_n_node(neighbour));  // U j-1
                                        }
                                    break;
                                case(2):
                                      if(bcs.get_bc(neighbour)){
                                            temp = 2*local_soln.get_rho(i)      // U j
                                                -soln.get_rho(neighbour)      // U j+1
                                                - soln.get_rho(Mesh.get_w_node(neighbour));

                                        }else{
                                            temp = soln.get_rho(Mesh.get_e_node(neighbour) ) // U j+2
                                                - 3* soln.get_rho(neighbour)                //U j+1
                                                + 3* local_soln.get_rho(i)                  // U j
                                            - soln.get_rho(Mesh.get_w_node(neighbour));  // U j-1
                                        }
                                    break;
                                case(3):
                                      if(bcs.get_bc(neighbour)){
                                            temp = 2*local_soln.get_rho(i)      // U j
                                                -soln.get_rho(neighbour)      // U j+1
                                                - soln.get_rho(Mesh.get_s_node(neighbour));

                                        }else{
                                            temp = soln.get_rho(Mesh.get_n_node(neighbour) ) // U j+2
                                                - 3* soln.get_rho(neighbour)                //U j+1
                                                + 3* local_soln.get_rho(i)                  // U j
                                                - soln.get_rho(Mesh.get_s_node(neighbour));  // U j-1
                                        }
                                    break;

                            }
                            break;

                        case(2):

                             switch(j){
                                case(0):
                                    if(bcs.get_bc(neighbour)){
                                            temp = 2*local_soln.get_u(i)      // U j
                                                -soln.get_u(neighbour)      // U j+1
                                                - soln.get_u(Mesh.get_e_node(neighbour));

                                    }else{
                                        temp = soln.get_u(Mesh.get_w_node(neighbour) ) // U j+2
                                        - 3* soln.get_u(neighbour)                //U j+1
                                        + 3* local_soln.get_u(i)                  // U j
                                        - soln.get_u(Mesh.get_e_node(neighbour));  // U j-1

                                    }

                                    break;
                                case(1):
                                      if(bcs.get_bc(neighbour)){
                                            temp = 2*local_soln.get_u(i)      // U j
                                                -soln.get_u(neighbour)      // U j+1
                                                - soln.get_u(Mesh.get_n_node(neighbour));

                                        }else{

                                        temp = soln.get_u(Mesh.get_s_node(neighbour) ) // U j+2
                                            - 3* soln.get_u(neighbour)                //U j+1
                                            + 3* local_soln.get_u(i)                  // U j
                                            - soln.get_u(Mesh.get_n_node(neighbour));  // U j-1
                                        }
                                    break;
                                case(2):
                                      if(bcs.get_bc(neighbour)){
                                            temp = 2*local_soln.get_u(i)      // U j
                                                -soln.get_u(neighbour)      // U j+1
                                                - soln.get_u(Mesh.get_w_node(neighbour));

                                        }else{
                                            temp = soln.get_u(Mesh.get_e_node(neighbour) ) // U j+2
                                                - 3* soln.get_u(neighbour)                //U j+1
                                                + 3* local_soln.get_u(i)                  // U j
                                            - soln.get_u(Mesh.get_w_node(neighbour));  // U j-1
                                        }
                                    break;
                                case(3):
                                      if(bcs.get_bc(neighbour)){
                                            temp = 2*local_soln.get_u(i)      // U j
                                                -soln.get_u(neighbour)      // U j+1
                                                - soln.get_u(Mesh.get_s_node(neighbour));

                                        }else{
                                            temp = soln.get_u(Mesh.get_n_node(neighbour) ) // U j+2
                                                - 3* soln.get_u(neighbour)                //U j+1
                                                + 3* local_soln.get_u(i)                  // U j
                                                - soln.get_u(Mesh.get_s_node(neighbour));  // U j-1
                                        }
                                    break;

                            }
                            break;
                        case(3):
                            switch(j){
                                case(0):
                                    if(bcs.get_bc(neighbour)){
                                            temp = 2*local_soln.get_v(i)      // U j
                                                -soln.get_v(neighbour)      // U j+1
                                                - soln.get_v(Mesh.get_e_node(neighbour));

                                    }else{
                                        temp = soln.get_v(Mesh.get_w_node(neighbour) ) // U j+2
                                        - 3* soln.get_v(neighbour)                //U j+1
                                        + 3* local_soln.get_v(i)                  // U j
                                        - soln.get_v(Mesh.get_e_node(neighbour));  // U j-1

                                    }

                                    break;
                                case(1):
                                      if(bcs.get_bc(neighbour)){
                                            temp = 2*local_soln.get_v(i)      // U j
                                                -soln.get_v(neighbour)      // U j+1
                                                - soln.get_v(Mesh.get_n_node(neighbour));

                                        }else{

                                        temp = soln.get_v(Mesh.get_s_node(neighbour) ) // U j+2
                                            - 3* soln.get_v(neighbour)                //U j+1
                                            + 3* local_soln.get_v(i)                  // U j
                                            - soln.get_v(Mesh.get_n_node(neighbour));  // U j-1
                                        }
                                    break;
                                case(2):
                                      if(bcs.get_bc(neighbour)){
                                            temp = 2*local_soln.get_v(i)      // U j
                                                -soln.get_v(neighbour)      // U j+1
                                                - soln.get_v(Mesh.get_w_node(neighbour));

                                        }else{
                                            temp = soln.get_v(Mesh.get_e_node(neighbour) ) // U j+2
                                                - 3* soln.get_v(neighbour)                //U j+1
                                                + 3* local_soln.get_v(i)                  // U j
                                            - soln.get_v(Mesh.get_w_node(neighbour));  // U j-1
                                        }
                                    break;
                                case(3):
                                      if(bcs.get_bc(neighbour)){
                                            temp = 2*local_soln.get_v(i)      // U j
                                                -soln.get_v(neighbour)      // U j+1
                                                - soln.get_v(Mesh.get_s_node(neighbour));

                                        }else{
                                            temp = soln.get_v(Mesh.get_n_node(neighbour) ) // U j+2
                                                - 3* soln.get_v(neighbour)                //U j+1
                                                + 3* local_soln.get_v(i)                  // U j
                                                - soln.get_v(Mesh.get_s_node(neighbour));  // U j-1
                                        }
                                    break;

                            }
                            break;
                        case(4):
                            switch(j){
                                case(0):
                                    if(bcs.get_bc(neighbour)){
                                            temp = 2*local_soln.get_w(i)      // U j
                                                -soln.get_w(neighbour)      // U j+1
                                                - soln.get_w(Mesh.get_e_node(neighbour));

                                    }else{
                                        temp = soln.get_w(Mesh.get_w_node(neighbour) ) // U j+2
                                        - 3* soln.get_w(neighbour)                //U j+1
                                        + 3* local_soln.get_w(i)                  // U j
                                        - soln.get_w(Mesh.get_e_node(neighbour));  // U j-1

                                    }

                                    break;
                                case(1):
                                      if(bcs.get_bc(neighbour)){
                                            temp = 2*local_soln.get_w(i)      // U j
                                                -soln.get_w(neighbour)      // U j+1
                                                - soln.get_w(Mesh.get_n_node(neighbour));

                                        }else{

                                        temp = soln.get_w(Mesh.get_s_node(neighbour) ) // U j+2
                                            - 3* soln.get_w(neighbour)                //U j+1
                                            + 3* local_soln.get_w(i)                  // U j
                                            - soln.get_w(Mesh.get_n_node(neighbour));  // U j-1
                                        }
                                    break;
                                case(2):
                                      if(bcs.get_bc(neighbour)){
                                            temp = 2*local_soln.get_w(i)      // U j
                                                -soln.get_w(neighbour)      // U j+1
                                                - soln.get_w(Mesh.get_w_node(neighbour));

                                        }else{
                                            temp = soln.get_w(Mesh.get_e_node(neighbour) ) // U j+2
                                                - 3* soln.get_w(neighbour)                //U j+1
                                                + 3* local_soln.get_w(i)                  // U j
                                            - soln.get_w(Mesh.get_w_node(neighbour));  // U j-1
                                        }
                                    break;
                                case(3):
                                      if(bcs.get_bc(neighbour)){
                                            temp = 2*local_soln.get_w(i)      // U j
                                                -soln.get_w(neighbour)      // U j+1
                                                - soln.get_w(Mesh.get_s_node(neighbour));

                                        }else{
                                            temp = soln.get_w(Mesh.get_n_node(neighbour) ) // U j+2
                                                - 3* soln.get_w(neighbour)                //U j+1
                                                + 3* local_soln.get_w(i)                  // U j
                                                - soln.get_w(Mesh.get_s_node(neighbour));  // U j-1
                                        }
                                    break;

                            }
                            break;


                    }

               return temp;

        }
