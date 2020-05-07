#include "Mesh.h"
#include <math.h>
#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include "global_variables.h"
#include <cstdio>
#include "TECIO.h"

#include <cstring>
#include <limits>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>




using namespace std;

Mesh::Mesh(domain_geometry domain, global_variables &globals)
{

    if( globals.mesh_type >2){

        return;
    }
    //ctor
    X = domain.X;
    Y = domain.Y;
    dx = domain.dx; //dimensional form
    dy = domain.dy;
    multi_grid_dt = domain.dt;

    if (globals.testcase == 3){
        num_x_cells = ceil(X/dx) + 4;
        num_y_cells= ceil(Y/dy) + 4;
    }else{
        num_x_cells = ceil(X/dx) + 2;
        num_y_cells= ceil(Y/dy) + 2;
	}
	// plus twos account for ghost nodes
	total_cells  = (num_x_cells ) * (num_y_cells);

    /// need error check here to see if grid divisible by multigrid criteria
    //Uniform only

	//dx = X/(num_x_cells-2); // reset dx/dy t0 allow for ceiling
    //dy = Y/(num_y_cells -2);
       // dx =dx/domain.dt; // turn into non-dimensional form
       // dy =dy/domain.dt;
        // 2.0 as cell is 2 LBM nodes in length
    cs = domain.cs;


    centroid_x = new double [total_cells +1];
        if (centroid_x==NULL) exit (1);
    centroid_y = new double [total_cells +1];
        if (centroid_y==NULL) exit (1);
    centroid_z = new double [total_cells +1];
        if (centroid_z==NULL) exit (1);
    north_x = new double [total_cells +1];
        if (north_x==NULL) exit (1);
    north_y = new double [total_cells +1];
        if (north_y==NULL) exit (1);
    north_z = new double [total_cells +1];
        if (north_z==NULL) exit (1);
    east_x = new double [total_cells +1];
        if (east_x==NULL) exit (1);
    east_y = new double [total_cells +1];
        if (east_y==NULL) exit (1);
    east_z = new double [total_cells +1];
        if (east_z==NULL) exit (1);
    west_x = new double [total_cells +1];
        if (west_x==NULL) exit (1);
    west_y = new double [total_cells +1];
        if (west_y==NULL) exit (1);
    west_z = new double [total_cells +1];
        if (west_z==NULL) exit (1);
    south_x = new double [total_cells +1];
        if (south_x==NULL) exit (1);
    south_y = new double [total_cells +1];
        if (south_y==NULL) exit (1);
    south_z = new double [total_cells +1];
        if (south_z==NULL) exit (1);

    n_area = new double [total_cells +1];
        if (n_area==NULL) exit (1);
    s_area = new double [total_cells +1];
        if (s_area==NULL) exit (1);
    w_area = new double [total_cells +1];
        if (w_area==NULL) exit (1);
    e_area = new double [total_cells +1];
        if (e_area==NULL) exit (1);

    cell_volume = new double [total_cells +1];
        if (cell_volume==NULL) exit (1);


    n_i = new double [total_cells +1];
        if (n_i==NULL) exit (1);
    n_j = new double [total_cells +1];
        if (n_j==NULL) exit (1);
    n_k = new double [total_cells +1];
        if (n_k==NULL) exit (1);
     e_i = new double [total_cells +1];
        if (e_i==NULL) exit (1);
    e_j = new double [total_cells +1];
        if (e_j==NULL) exit (1);
    e_k = new double [total_cells +1];
        if (e_k==NULL) exit (1);
     w_i = new double [total_cells +1];
        if (w_i==NULL) exit (1);
    w_j = new double [total_cells +1];
        if (w_j==NULL) exit (1);
    w_k = new double [total_cells +1];
        if (w_k==NULL) exit (1);
     s_i = new double [total_cells +1];
        if (s_i==NULL) exit (1);
    s_j = new double [total_cells +1];
        if (s_j==NULL) exit (1);
    s_k = new double [total_cells +1];
        if (s_k==NULL) exit (1);


     n_node = new int [total_cells +1];
        if (n_node==NULL) exit (1);
    s_node = new int [total_cells +1];
        if (s_node==NULL) exit (1);
    w_node = new int [total_cells +1];
        if (w_node==NULL) exit (1);
    e_node = new int [total_cells +1];
        if (e_node==NULL) exit (1);




    delta_t = new double [total_cells +1];
        if (delta_t==NULL) exit (1);
    delta_t_n = new double [total_cells +1];
    if (delta_t==NULL) exit (1);
    delta_t_e = new double [total_cells +1];
    if (delta_t==NULL) exit (1);

    if (globals.mesh_type == 1) {
        if(globals.testcase == 3){
            this->create_standard_mesh(-1,-1);
        }else{
            this->create_standard_mesh(0,0);
        }


    }else if(globals.mesh_type == 2){
        this -> create_cosine_mesh(globals.PI, domain.dt);
    }



        // allocate memory to arrays

    std::ofstream mesh_output ;
    std::string output_dir;

    output_dir = globals.output_file +"/mesh.txt";


        // error_output.open("/home/brendan/Dropbox/PhD/Test Cases/Couette Flow/error.txt", ios::out);
    mesh_output.open(output_dir.c_str(), ios::out);
    for (int mesh_t= 0; mesh_t < total_cells; mesh_t++){
         mesh_output << mesh_t << "," << centroid_x[mesh_t] << "," <<
                centroid_y[mesh_t] << "," << centroid_z[mesh_t]
                << ","
                << delta_t[mesh_t]  << ","
                << delta_t_n[mesh_t] << ","<< delta_t_e[mesh_t]
                << "," <<  (globals.visc *3/delta_t[mesh_t]) << endl;
        }

    mesh_output.close();

    globals.update_visc(Y);


	X=X;
}

Mesh::~Mesh()
{
    //dtor
    delete [](centroid_x);
    centroid_x = NULL;
    delete [] (north_x);
    north_x = NULL;
    delete [] (south_x);
    south_x = NULL;
    delete [] (east_x);
    east_x = NULL;
    delete [] (west_x);
    west_x = NULL;


    delete [](centroid_y);
    centroid_y = NULL;
    delete [] (north_y);
    north_y = NULL;
    delete [] (south_y);
    south_y = NULL;
    delete [] (east_y);
    east_y = NULL;
    delete [] (west_y);
    west_y = NULL;

    delete [](centroid_z);
    centroid_z = NULL;
    delete [] (north_z);
    north_z = NULL;
    delete [] (south_z);
    south_z = NULL;
    delete [] (east_z);
    east_z = NULL;
    delete [] (west_z);
    west_z = NULL;

    delete [](cell_volume);
    cell_volume = NULL;
    delete [](delta_t);
    delta_t = NULL;
    delete [](delta_t_n);
    delta_t_n = NULL;
    delete [](delta_t_e);
    delta_t_e = NULL;

    delete [](n_area);
    n_area = NULL;
    delete [](n_i);
    n_i = NULL;
    delete [](n_j);
    n_j = NULL;
    delete [](n_k);
    n_k = NULL;
    delete [](n_node);
    n_node = NULL;

    delete [](e_area);
    e_area = NULL;
    delete [](e_i);
    e_i = NULL;
    delete [](e_j);
    e_j = NULL;
    delete [](e_k);
    e_k = NULL;
    delete [](e_node);
    e_node = NULL;

    delete [](w_area);
    w_area = NULL;
    delete [](w_i);
    w_i = NULL;
    delete [](w_j);
    w_j = NULL;
    delete [](w_k);
    w_k = NULL;
    delete [](w_node);
    w_node = NULL;

    delete [](s_area);
    s_area = NULL;
    delete [](s_i);
    s_i = NULL;
    delete [](s_j);
    s_j = NULL;
    delete [](s_k);
    s_k = NULL;
    delete [](s_node);
    s_node = NULL;
}

void Mesh::create_standard_mesh(int strt, int fnsh ){
    int counter =0;
    for( int j=(0 +strt); j < (num_y_cells + fnsh); j++){
        for( int i=(0+strt); i < (num_x_cells + fnsh); i++){

            centroid_x[counter] = dx/2 + (i-1)*dx;
            centroid_y[counter] = dy/2 + (j-1)*dy;
            centroid_z[counter] = 0; //temporary
            north_x[counter] = dx/2 + (i-1)*dx;
            north_y[counter] = dy + (j-1)*dy;
            north_z[counter] = 0; //temporary
            south_x[counter] = dx/2 + (i-1)*dx;
            south_y[counter] = (j-1)*dy;
            south_z[counter] = 0; //temporary
            west_x[counter] = (i-1)*dx;
            west_y[counter] = dy/2 + (j-1)*dy;
            west_z[counter] = 0; //temporary
            east_x[counter] = dx + (i-1)*dx;
            east_y[counter] = dy/2 + (j-1)*dy;
            east_z[counter] = 0; //temporary

            n_area[counter] = dx;
            s_area[counter] = dx;
            e_area[counter] = dy;
            w_area[counter] = dy;

            cell_volume[counter] = dx*dy;

            n_i[counter] = 0.0;
            n_j[counter] = 1.0;
            n_k[counter] = 0; //temporary

            e_i[counter] = 1.0;
            e_j[counter] = 0.0;
            e_k[counter] = 0; //temporary


            s_i[counter] = 0.0;
            s_j[counter] = -1.0;
            s_k[counter] = 0; //temporary


            w_i[counter] = -1.0;
            w_j[counter] = 0.0;
            w_k[counter] = 0; //temporary

            delta_t[counter] = 0.5 * std::min(dy,dx);
            delta_t_n[counter] = delta_t[counter];
            delta_t_e[counter] = delta_t[counter];

             // West boundary
            if( i ==(0 + strt)){
                w_node[counter] = -1 ;//dummy value
            }else{
                //w_node[counter] = counter - num_y_cells;
                w_node[counter] = counter - 1;
            }

            // east boundary
            if ( i == (num_x_cells -1 + fnsh)){
                e_node[counter] = -1; //dummy value
            }else{
                e_node[counter] = counter + 1;
            }

            // south boundary
            if(j == (0+ strt) ){
                s_node[counter] = -1; //dummy value
            }else{
                s_node[counter] = counter - num_x_cells ;
            }

            // north boundary
            if( j == (num_y_cells-1 +fnsh)){
                n_node[counter] = -1; //dummy value
            }else{
                n_node[counter] = counter + num_x_cells ;
            }

            counter = counter +1;
        }
	}
    int i;
	// update corner cells neighbours to not include

	//SW
	i =1;
	w_node[i] = -1;
	i = num_x_cells;
	s_node[i] = -1;

	//SE
	i = num_x_cells -2;
	e_node[i] = -1;
	i = 2 *num_x_cells -1;
	s_node[i] = -1;


	//NW
	i = (num_y_cells - 2) * num_x_cells ;
    n_node[i] = -1;
    i = (num_y_cells - 1) * num_x_cells +1;
	w_node[i] = -1;

	//NE
	i = (num_y_cells -1 ) * num_x_cells - 1;
    n_node[i] = -1;
    i = (num_y_cells ) * num_x_cells - 2;
	e_node[i] = -1;

}





void Mesh::calc_unstructured_lattice_size(){


    for (int i = 0; i< total_cells; i++){
            //get the north and east flux streaming time step
            delta_t_n[i] = min(delta_t[i], delta_t[n_node[i]]);
            delta_t_e[i] = min(delta_t[i], delta_t[e_node[i]]);
    }

}

void Mesh::get_centroid(int i, vector_var &cell){

    cell.x = centroid_x[i];
    cell.y = centroid_y[i];
    cell.z = centroid_z[i];

}

void Mesh::create_cosine_mesh(double PI,double dt){
    int counter =0;

    double delta_y ,delta_x , y2,y1,x2,x1,cx,cy;
    double min_xy;
    min_xy = 1/(0.5* ( 1- cos((1)/(Y) * PI)));

    for( int j=0; j < num_y_cells ; j++){
        for( int i=0; i < num_x_cells; i++){
            y2 = 0.5* ( 1- cos((j+1-1)/(Y ) * PI))*min_xy;
            y1 = 0.5* ( 1- cos((j-1)/(Y ) * PI))*min_xy;
            x2 = 0.5* ( 1- cos((i+1-1)/(X ) * PI))*min_xy;
            x1 = 0.5* ( 1- cos((i-1)/(X ) * PI))*min_xy;
            if(j ==0){
                y1 = y1*-1;

            }
            if(i ==0){
                x1 = x1*-1;
            }
            if(j == num_y_cells-1){
                y2 = min_xy + min_xy- y2;
            }
            if(i == num_x_cells -1){
                x2  = min_xy + min_xy -x2;
            }

            delta_y = (y2-y1);
            delta_x = (x2-x1);
            cx = (x1+x2)/2;
            cy = (y1+y2)/2;

            centroid_x[counter] = cx;
            centroid_y[counter] = cy;
            centroid_z[counter] = 0; //temporary
            north_x[counter] = cx ;
            north_y[counter] = y2;
            north_z[counter] = 0; //temporary
            south_x[counter] = cx;
            south_y[counter] = y1;
            south_z[counter] = 0; //temporary
            west_x[counter] = x1;
            west_y[counter] = cy;
            west_z[counter] = 0; //temporary
            east_x[counter] = x2;
            east_y[counter] = cy;
            east_z[counter] = 0; //temporary

            n_area[counter] = delta_x;
            s_area[counter] = delta_x;
            e_area[counter] = delta_y;
            w_area[counter] = delta_y;

            cell_volume[counter] = delta_x*delta_y;

            n_i[counter] = 0.0;
            n_j[counter] = 1.0;
            n_k[counter] = 0; //temporary

            e_i[counter] = 1.0;
            e_j[counter] = 0.0;
            e_k[counter] = 0; //temporary


            s_i[counter] = 0.0;
            s_j[counter] = -1.0;
            s_k[counter] = 0; //temporary


            w_i[counter] = -1.0;
            w_j[counter] = 0.0;
            w_k[counter] = 0; //temporary

            delta_t[counter] = 0.5 * std::min(delta_x,delta_y);

             // West boundary
            if( i ==0){
                w_node[counter] = -1 ;//dummy value
            }else{
                //w_node[counter] = counter - num_y_cells;
                w_node[counter] = counter - 1;
            }

            // east boundary
            if ( i == (num_x_cells -1)){
                e_node[counter] = -1; //dummy value
            }else{
                e_node[counter] = counter + 1;
            }

            // south boundary
            if(j == 0){
                s_node[counter] = -1; //dummy value
            }else{
                s_node[counter] = counter - num_x_cells;
            }

            // north boundary
            if( j == (num_y_cells-1)){
                n_node[counter] = -1; //dummy value
            }else{
                n_node[counter] = counter + num_x_cells;
            }

            counter = counter +1;
        }
	}

    int i;
	// update corner cells
	//SW
	i =1;
	w_node[i] = -1;
	i = num_x_cells;
	s_node[i] = -1;

	//SE
	i = num_x_cells -2;
	e_node[i] = -1;
	i = 2 *num_x_cells -1;
	s_node[i] = -1;


	//NW
	i = (num_y_cells - 2) * num_x_cells ;
    n_node[i] = -1;
    i = (num_y_cells - 1) * num_x_cells +1;
	w_node[i] = -1;

	//NE
	i = (num_y_cells -1 ) * num_x_cells - 1;
    n_node[i] = -1;
    i = (num_y_cells ) * num_x_cells - 2;
	e_node[i] = -1;

	  X = min_xy;
    Y = min_xy;
    calc_unstructured_lattice_size();

}


double  Mesh::get_node_x(int node_num){
    double result;
    result = centroid_x[node_num];

    return result ;
}



domain_geometry Mesh::create_coarse_mesh_domain(){

    domain_geometry coarse_domain;

    coarse_domain.X = X/dx;
    coarse_domain.Y = Y/dy;
    coarse_domain.dx = dx*2*multi_grid_dt;
    coarse_domain.dy = dy*2*multi_grid_dt;
    coarse_domain.dt = multi_grid_dt; //streaming time step
    //coarse_domain.dt = multi_grid_dt;
    coarse_domain.cs = cs;


    return coarse_domain;


}




