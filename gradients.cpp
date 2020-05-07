#include "gradients.h"


gradients::gradients()
{
    //ctor
}


gradients::gradients(int _total_nodes)
{
    //ctor
    total_nodes = _total_nodes;
     //rho = (double*) malloc (sizeof(double)*(total_nodes));
     u = new vector_var [total_nodes +1];
        if (u==NULL) exit (1);
     v = new vector_var [total_nodes +1];
        if (v==NULL) exit (1);
     w = new vector_var [total_nodes +1];
        if (w==NULL) exit (1);
     rho = new vector_var [total_nodes +1];
        if (rho==NULL) exit (1);
		
		LHS_xx = new double[total_nodes + 1];
		if (LHS_xx == NULL) exit(1);
		
		LHS_xy = new double[total_nodes + 1];
		if (LHS_xy == NULL) exit(1);
		
		LHS_xz = new double[total_nodes + 1];
		if (LHS_xz == NULL) exit(1);

		LHS_yx = new double[total_nodes + 1];
		if (LHS_yx == NULL) exit(1);

		LHS_yy = new double[total_nodes + 1];
		if (LHS_yy == NULL) exit(1);

		LHS_yz = new double[total_nodes + 1];
		if (LHS_yz == NULL) exit(1);

		LHS_zx = new double[total_nodes + 1];
		if (LHS_zx == NULL) exit(1);

		LHS_zy = new double[total_nodes + 1];
		if (LHS_zy == NULL) exit(1);

		LHS_zz = new double[total_nodes + 1];
		if (LHS_zz == NULL) exit(1);

		RHS_x = new double[6*(total_nodes + 1)];
		if (RHS_x == NULL) exit(1);

		RHS_y = new double[6 * (total_nodes + 1)];
		if (RHS_y == NULL) exit(1);

		RHS_z = new double[6 * (total_nodes + 1)];
		if (RHS_z == NULL) exit(1);
		

    Initialise();

}

gradients::~gradients()
{
    //dtor
    delete [] u;
    u = NULL;
    delete [] v;
    v = NULL;
    delete [] w;
    w = NULL;
    delete [] rho;
    rho = NULL;

	delete[] LHS_xx;
	LHS_xx = NULL;
	delete[] LHS_xy;
	LHS_xy = NULL;
	delete[] LHS_xz;
	LHS_xz = NULL;

	delete[] LHS_yx;
	LHS_yx = NULL;
	delete[] LHS_yy;
	LHS_yy = NULL;
	delete[] LHS_yz;
	LHS_yz = NULL;

	delete[] LHS_zx;
	LHS_zx = NULL;
	delete[] LHS_zy;
	LHS_zy = NULL;
	delete[] LHS_zz;
	LHS_zz = NULL;

	delete[] RHS_x;
	RHS_x = NULL;
	delete[] RHS_y;
	RHS_y = NULL;
	delete[] RHS_z;
	RHS_z = NULL;

}

void gradients::Initialise() {

   for( int i =0; i< total_nodes; i++){
        u[i].zero();
        v[i].zero();
        w[i].zero();
        rho[i].zero();

    }
}


void gradients::add_LS_contributions(vector_var &RHS_rho, vector_var &RHS_u, vector_var &RHS_v, vector_var &RHS_w, Boundary_Conditions &bcs,unstructured_mesh &mesh,domain_geometry &domain,
                      Solution &src , int i1, int i,int i_nb) {

        double d_u,d_v,d_w,d_rho;

       //boundary condition, find gradient at shared face

        d_u = (src.get_u(i1) - src.get_u(i));
        d_v = (src.get_v(i1) - src.get_v(i));
        d_w = (src.get_w(i1) - src.get_w(i));
        d_rho = (src.get_rho(i1) - src.get_rho(i)) ;

		RHS_rho.x = RHS_rho.x + RHS_x[i * 6 +i_nb] * d_rho;
		RHS_rho.y = RHS_rho.y + RHS_y[i * 6 + i_nb] * d_rho;
		RHS_rho.z = RHS_rho.z + RHS_z[i * 6 + i_nb] * d_rho;

		RHS_u.x = RHS_u.x + RHS_x[i * 6 + i_nb] * d_u;
		RHS_u.y = RHS_u.y + RHS_y[i * 6 + i_nb] * d_u;
		RHS_u.z = RHS_u.z + RHS_z[i * 6 + i_nb] * d_u;

		RHS_v.x = RHS_v.x + RHS_x[i * 6 + i_nb] * d_v;
		RHS_v.y = RHS_v.y + RHS_y[i * 6 + i_nb] * d_v;
		RHS_v.z = RHS_v.z + RHS_z[i * 6 + i_nb] * d_v;

		RHS_w.x = RHS_w.x + RHS_x[i * 6 + i_nb] * d_w;
		RHS_w.y = RHS_w.y + RHS_y[i * 6 + i_nb] * d_w;
		RHS_w.z = RHS_w.z + RHS_z[i * 6 + i_nb] * d_w;

}

void gradients::add_non_macrovariable_contributions(Boundary_Conditions &bcs,
                        unstructured_mesh &mesh,domain_geometry &domain,
                      Solution &src , int i1, int i,double beta_0, int face,int i_nb) {

        double dx,dy,dz,w;
        double alpha_w;
        double L,L_norm, l_norm,s;
        


         // distance along fromcell centre to face centre

        dx = mesh.get_face_x(face) - mesh.get_centroid_x(i);
        dy = mesh.get_face_y(face) - mesh.get_centroid_y(i);
        dz = mesh.get_face_z(face) - mesh.get_centroid_z(i);

         // distance along normal to face centre
        l_norm = dx* mesh.get_face_i(face) + dy * mesh.get_face_j(face) + dz * mesh.get_face_k(face);

        //get distance to neighbouring cell centre
        dx = mesh.get_centroid_x(i1) - mesh.get_centroid_x(i);
        dy = mesh.get_centroid_y(i1) - mesh.get_centroid_y(i);
        dz = mesh.get_centroid_z(i1) - mesh.get_centroid_z(i);

        // L distance between two cell centres
        L = sqrt(pow(dx,2.0) + pow(dy,2.0) +pow(dz,2.0) );

        // area of interface
        s = mesh.get_face_area(face);

        // distance along normal to cell centre
        L_norm = dx* mesh.get_face_i(face) + dy * mesh.get_face_j(face) + dz * mesh.get_face_k(face);

        alpha_w = pow( 2* l_norm/L_norm ,2.0);

        w = alpha_w * s/L;
		       
		LHS_xx[i] = LHS_xx[i] + w * dx*dx*beta_0;
		LHS_xy[i] = LHS_xy[i] + w * dx*dy*beta_0;
		LHS_xz[i] = LHS_xz[i] + w * dx*dz*beta_0;
		LHS_yx[i] = LHS_yx[i] + w * dy*dx*beta_0;
		LHS_yy[i] = LHS_yy[i] + w * dy*dy*beta_0;
		LHS_yz[i] = LHS_yz[i] + w * dy*dz*beta_0;
		LHS_zx[i] = LHS_zx[i] + w * dz*dx*beta_0;
		LHS_zy[i] = LHS_zy[i] + w * dz*dy*beta_0;
		LHS_zz[i] = LHS_zz[i] + w * dz*dz*beta_0;



        //assumed that alpha G,j is 0.5 to give true Green gauss formula when beta_0 =0
        if(mesh.get_mesh_owner(face) == i){
            RHS_x[i*6+ i_nb] = (beta_0 * w *dx + 2*(1- beta_0)*0.5  *mesh.get_face_i(face)  * mesh.get_face_area(face));
			RHS_y[i * 6 + i_nb] =  (beta_0 * w *dy + 2*(1- beta_0)*0.5  *mesh.get_face_j(face) * mesh.get_face_area(face)) ;
			RHS_z[i * 6 + i_nb] = (beta_0 * w *dz + 2*(1- beta_0)*0.5  *mesh.get_face_k(face) * mesh.get_face_area(face)) ;
        }else{
			RHS_x[i * 6 + i_nb] = (beta_0 * w *dx + 2*(1- beta_0)*0.5  *mesh.get_face_i(face) *-1 * mesh.get_face_area(face));
			RHS_y[i * 6 + i_nb] =  (beta_0 * w *dy + 2*(1- beta_0)*0.5  *mesh.get_face_j(face) *-1 * mesh.get_face_area(face)) ;
			RHS_z[i * 6 + i_nb] = (beta_0 * w *dz + 2*(1- beta_0)*0.5  *mesh.get_face_k(face) *-1 * mesh.get_face_area(face)) ;

        }
		w = w;
}


//update gradients for each cell
void gradients::pre_fill_LHS_and_RHS_matrix(Boundary_Conditions &bcs,unstructured_mesh &mesh,domain_geometry &domain,
                      Solution &src, global_variables &globals ){

   double beta_0,dx,dy,dz;
    int face,i1;

    double max_disp,max_area,disp, aspect_ratio;
    max_area = 0.0;
    max_disp =0.0;
	double a, b, c, d, e, f, g, h, k;
	double det_LHS;
   for(int i = 0; i < mesh.get_n_cells(); i++){
	   det_LHS = 0.0;
	   LHS_xx[i] = 0.0;
	   LHS_xy[i] = 0.0;
	   LHS_xz[i] = 0.0;
	   LHS_yx[i] = 0.0;
	   LHS_yy[i] = 0.0;
	   LHS_yz[i] = 0.0;
	   LHS_zx[i] = 0.0;
	   LHS_zy[i] = 0.0;
	   LHS_zz[i] = 0.0;
	   max_area = 0.0;
	   max_disp = 0.0;

	   //loop through all faces to find max area and displacement
         for (int i_nb =0 ; i_nb < mesh.gradient_cells[i].size(); i_nb++){
            i1 = mesh.gradient_cells[i][i_nb];
            face = mesh.gradient_faces[i][i_nb];

            max_area = std::max(max_area, mesh.get_face_area(face));

            dx = mesh.get_centroid_x(i1) - mesh.get_centroid_x(i);
            dy = mesh.get_centroid_y(i1) - mesh.get_centroid_y(i);
            dz = mesh.get_centroid_z(i1) - mesh.get_centroid_z(i);

            disp = sqrt(dx*dx + dy*dy + dz*dz );
            max_disp = std::max(disp,max_disp);
        }

        //from Shima
        aspect_ratio = 2 * fabs(max_disp) * max_area/mesh.get_cell_volume(i);
       
        //get beta factor for each cell
		//beta is tuner for how GG or LS the hyrbid approach should be
        if(globals.gradient_calc_type == "hybrid"){
            beta_0 = std::min(1.0, (2 / aspect_ratio));
        }else if(globals.gradient_calc_type == "LS"){
            beta_0 = 1.0;
        }else{ // green gauss
            beta_0 = 0.0;
        }
       
		LHS_xx[i] = LHS_xx[i] + 2*(1-beta_0) *mesh.get_cell_volume(i);
		LHS_yy[i] = LHS_yy[i] + 2*(1-beta_0) *mesh.get_cell_volume(i);
		LHS_zz[i] = LHS_zz[i] + 2*(1-beta_0) *mesh.get_cell_volume(i);

        //calculate LHS and RHS contributions
        for (int i_nb =0 ; i_nb < mesh.gradient_cells[i].size(); i_nb++){
            i1 = mesh.gradient_cells[i][i_nb];
            face = mesh.gradient_faces[i][i_nb];

            if(i1 < 0){
                i1 = i1;
            }
            add_non_macrovariable_contributions(bcs,mesh,domain,src,i1,i,beta_0,face, i_nb);

        }
		det_LHS = LHS_xx[i] * (LHS_yy[i] * LHS_zz[i] - LHS_yz[i] * LHS_zy[i]) +
			LHS_xy[i] * (LHS_yx[i] * LHS_zz[i] - LHS_yz[i] - LHS_zx[i]) +
			LHS_xz[i]* (LHS_yx[i] * LHS_zy[i] - LHS_yy[i] * LHS_zx[i]);


		//temp variables
		a = LHS_xx[i];
		b = LHS_xy[i];
		c = LHS_xz[i];
		d = LHS_yx[i];
		e = LHS_yy[i];
		f = LHS_yz[i];
		g = LHS_zx[i];
		h = LHS_zy[i];
		k = LHS_zz[i];

		//store inverse
		LHS_xx[i] = 1 / det_LHS * (e*k - f * h);
		LHS_xy[i] = 1 / det_LHS * (c*h -b*k);
		LHS_xz[i] = 1 / det_LHS * (b*f-c*e);
		LHS_yx[i] = 1 / det_LHS * (f*g-d*k);
		LHS_yy[i] = 1 / det_LHS * (a*k-c*g);
		LHS_yz[i] = 1 / det_LHS * (c*d-a*f);
		LHS_zx[i] = 1 / det_LHS * (d*h-e*g);
		LHS_zy[i] = 1 / det_LHS * (b*g-a*h);
		LHS_zz[i] = 1 / det_LHS * (a*e-b*d);

   }



  }




//update gradients for each cell
void gradients::Get_LS_Gradients(Boundary_Conditions &bcs,unstructured_mesh &mesh,domain_geometry &domain,
                      Solution &src, global_variables &globals ){

    int i1 ;
    int neighbour;
    vector_var RHS_u,RHS_v,RHS_w,RHS_rho;
	vector_var grad_u, grad_v, grad_w, grad_rho;
    vector_var cell_to_face,grad_cell,bc_correction;
    //interior cells

    cell_to_face.set_equal(0.0,0.0,0.0);
    bc_correction.set_equal(0.0,0.0,0.0);
    grad_cell.set_equal(0.0,0.0,0.0);


//    double min_u,max_u, min_rho, max_rho;
//
//    max_u = 2*globals.max_velocity*0.1;
//    min_u = -max_u;
//
//    min_rho = -0.003;
//    max_rho = 0.003;
    for(int i =0; i< mesh.get_n_cells();i++){

        // reset Matrix


		RHS_u.set_equal(0.0, 0.0, 0.0);
        RHS_v.set_equal(0.0, 0.0, 0.0);
        RHS_w.set_equal(0.0, 0.0, 0.0);
        RHS_rho.set_equal(0.0, 0.0, 0.0);
        // delta distance

        for (int i_nb =0 ; i_nb < mesh.gradient_cells[i].size(); i_nb++){
            i1 = mesh.gradient_cells[i][i_nb];

            if(i1 < 0){
                i1 = i1;
            }
            add_LS_contributions(RHS_rho,RHS_u,RHS_v,RHS_w ,bcs,mesh,domain,src,i1,i,i_nb);

        }
		
		
		rho[i].x = LHS_xx[i] * RHS_rho.x + LHS_xy[i] * RHS_rho.y + LHS_xz[i] * RHS_rho.z;
		rho[i].y = LHS_yx[i] * RHS_rho.x + LHS_yy[i]  * RHS_rho.y + LHS_yz[i] * RHS_rho.z;
		rho[i].z = LHS_zx[i] * RHS_rho.x + LHS_zy[i] * RHS_rho.y + LHS_zz[i]  * RHS_rho.z;

		u[i].x = LHS_xx[i] * RHS_u.x + LHS_xy[i] * RHS_u.y + LHS_xz[i] * RHS_u.z;
		u[i].y = LHS_yx[i] * RHS_u.x + LHS_yy[i] * RHS_u.y + LHS_yz[i] * RHS_u.z;
		u[i].z = LHS_zx[i] * RHS_u.x + LHS_zy[i] * RHS_u.y + LHS_zz[i] * RHS_u.z;

		v[i].x = LHS_xx[i] * RHS_v.x + LHS_xy[i] * RHS_v.y + LHS_xz[i] * RHS_v.z;
		v[i].y = LHS_yx[i] * RHS_v.x + LHS_yy[i] * RHS_v.y + LHS_yz[i] * RHS_v.z;
		v[i].z = LHS_zx[i] * RHS_v.x + LHS_zy[i] * RHS_v.y + LHS_zz[i] * RHS_v.z;

		w[i].x = LHS_xx[i] * RHS_w.x + LHS_xy[i] * RHS_w.y + LHS_xz[i] * RHS_w.z;
		w[i].y = LHS_yx[i] * RHS_w.x + LHS_yy[i] * RHS_w.y + LHS_yz[i] * RHS_w.z;
		w[i].z = LHS_zx[i] * RHS_w.x + LHS_zy[i] * RHS_w.y + LHS_zz[i] * RHS_w.z;

		
        }

        int j;


    ///BCS will find grads on face centre of shared face


		RHS_u.set_equal(0.0, 0.0, 0.0);
		RHS_v.set_equal(0.0, 0.0, 0.0);
		RHS_w.set_equal(0.0, 0.0, 0.0);
		RHS_rho.set_equal(0.0, 0.0, 0.0);

    int face;

    for( int i =0; i < mesh.get_num_bc(); i++){
        cell_to_face.set_equal(0.0,0.0,0.0);

         face = mesh.get_n_neighbours() + i;
        neighbour = mesh.get_mesh_owner(face);
        //periodic nodes get LS treatment
        j = i + mesh.get_n_cells();


        /// assign gradients

        //dirichlet
        if(bcs.get_rho_type(i) == 1){
            get_ghost_grads(mesh, src,i,neighbour,true,bcs,face,globals);

        //neumann
        }else if (bcs.get_rho_type(i) == 2){
            rho[j].x = rho[neighbour].x *(1- mesh.get_face_i(face));
            rho[j].y = rho[neighbour].y*(1- mesh.get_face_j(face));
            rho[j].z = rho[neighbour].z*(1- mesh.get_face_k(face));

       //inlet conditions
        }else if (bcs.get_rho_type(i) == 6){

            rho[j].x = 0.0;
            rho[j].y = 0.0;
            rho[j].z = 0.0;

         //wall
        }else if (bcs.get_rho_type(i) == 7){
            rho[j].x = rho[neighbour].x ;
            rho[j].y = rho[neighbour].y;
            rho[j].z = rho[neighbour].z;
             //periodic are full LS

        }else if (bcs.get_rho_type(i) == 8){
             get_ghost_grads(mesh, src,i,neighbour,false,bcs,face,globals);
             //periodic are full LS

        }else {

            rho[j].x = rho[bcs.get_periodic_node(i)].x;
            rho[j].y = rho[bcs.get_periodic_node(i)].y;
            rho[j].z = rho[bcs.get_periodic_node(i)].z;

        }


         //dirichlet
        if(bcs.get_vel_type(i) == 1){
            get_ghost_grads(mesh, src,i,neighbour,false,bcs,face,globals);

        //neumann
        }else if (bcs.get_vel_type(i) == 2){

            // page 612 of Moukallad
            u[j].x = u[neighbour].x *(1- mesh.get_face_i(face));
            u[j].y = u[neighbour].y*(1- mesh.get_face_j(face));
            u[j].z = u[neighbour].z*(1- mesh.get_face_k(face));

            v[j].x = v[neighbour].x *(1- mesh.get_face_i(face));
            v[j].y = v[neighbour].y*(1- mesh.get_face_j(face));
            v[j].z = v[neighbour].z*(1- mesh.get_face_k(face));


            w[j].x = w[neighbour].x *(1- mesh.get_face_i(face));
            w[j].y = w[neighbour].y*(1- mesh.get_face_j(face));
            w[j].z = w[neighbour].z*(1- mesh.get_face_k(face));

        //periodic are full LS

        //inlet conditions
        }else if (bcs.get_vel_type(i) == 6){

            u[j].x = 0.0;
            u[j].y =0.0;
            u[j].z = 0.0;

            v[j].x = 0.0;
            v[j].y = 0.0;
            v[j].z = 0.0;

            w[j].x = 0.0;
            w[j].y = 0.0;
            w[j].z = 0.0;
        }else if (bcs.get_vel_type(i) == 8){
            //symmetry about x axis
            get_ghost_grads(mesh, src,i,neighbour,false,bcs,face,globals);

        //periodic are full LS

        //inlet conditions
        }else {

            u[j].x = u[bcs.get_periodic_node(i)].x;
            u[j].y = u[bcs.get_periodic_node(i)].y;
            u[j].z = u[bcs.get_periodic_node(i)].z;

            v[j].x = v[bcs.get_periodic_node(i)].x;
            v[j].y = v[bcs.get_periodic_node(i)].y;
            v[j].z = v[bcs.get_periodic_node(i)].z;

            w[j].x = w[bcs.get_periodic_node(i)].x;
            w[j].y = w[bcs.get_periodic_node(i)].y;
            w[j].z = w[bcs.get_periodic_node(i)].z;


        }

    }

}
double gradients::limit_grad(double min_u, double max_u, double grad_u ){


    if(grad_u > 0){
        if( grad_u > max_u){
            return max_u;
        }else{
            return grad_u;
        }

    }else{
            if( grad_u < min_u){
                return min_u;
            }else{
                return grad_u;
            }

    }

}
void gradients::get_ghost_grads_limits( unstructured_mesh &mesh,Solution &src , int bc, int neighbour,bool rho_type,
                        Boundary_Conditions &bcs,int face,global_variables &globals) {

           double min_u,max_u, min_rho, max_rho;

        max_u = 2*globals.max_velocity;
        min_u = -max_u;

        min_rho = -0.003;
        max_rho = 0.003;

        double n_i,n_j,n_k;


       double dx,dy,dz,d_u,d_v,d_w,d_rho;
        int i1;
        i1 = bc + mesh.get_n_cells();
        int i_face;
        i_face = bc + mesh.get_n_neighbours();

        double d_mag;

       dx = mesh.get_centroid_x(i1)- mesh.get_centroid_x(neighbour);
       dy = mesh.get_centroid_y(i1) - mesh.get_centroid_y(neighbour);
       dz = mesh.get_centroid_z(i1) - mesh.get_centroid_z(neighbour);

        d_mag = sqrt(pow(dx,2.0) + pow(dy,2.0) + pow(dz,2.0));

        d_u = 2*(bcs.get_u(bc)- src.get_u(neighbour)) ;
        d_v = 2*(bcs.get_v(bc) - src.get_v(neighbour));
        d_w = 2*(bcs.get_w(bc) - src.get_w(neighbour));
        d_rho = 2*(bcs.get_rho(bc) - src.get_rho(neighbour)) ;

        n_i = mesh.get_face_i(i_face);
        n_j = mesh.get_face_j(i_face);
        n_k = mesh.get_face_k(i_face);

        if( rho_type){


            rho[i1].x = limit_grad(min_rho,max_rho,d_rho/d_mag*n_i);
            rho[i1].y = limit_grad(min_rho,max_rho,d_rho/d_mag*n_j);
            rho[i1].z = limit_grad(min_rho,max_rho,d_rho/d_mag*n_k);


        }else{
                u[i1].x = limit_grad(min_u,max_u,d_u/d_mag*n_i);
                v[i1].x = limit_grad(min_u,max_u,d_v/d_mag*n_i);
                w[i1].x = limit_grad(min_u,max_u,d_w/d_mag*n_i);

                u[i1].y = limit_grad(min_u,max_u,d_u/d_mag*n_j);
                v[i1].y = limit_grad(min_u,max_u,d_v/d_mag*n_j);
                w[i1].y = limit_grad(min_u,max_u,d_w/d_mag*n_j);

                u[i1].z = limit_grad(min_u,max_u,d_u/d_mag*n_k);
                v[i1].z = limit_grad(min_u,max_u,d_v/d_mag*n_k);
                w[i1].z = limit_grad(min_u,max_u,d_w/d_mag*n_k);
            }




    }


void gradients::get_ghost_grads( unstructured_mesh &mesh,Solution &src , int bc, int neighbour,bool rho_type,
                        Boundary_Conditions &bcs,int face,global_variables &globals) {



        double n_i,n_j,n_k;


       double dx,dy,dz,d_u,d_v,d_w,d_rho;
        int i1;
        i1 = bc + mesh.get_n_cells();
        int i_face;
        i_face = bc + mesh.get_n_neighbours();

        double d_mag;

       dx = mesh.get_centroid_x(i1)- mesh.get_centroid_x(neighbour);
       dy = mesh.get_centroid_y(i1) - mesh.get_centroid_y(neighbour);
       dz = mesh.get_centroid_z(i1) - mesh.get_centroid_z(neighbour);

        d_mag = sqrt(pow(dx,2.0) + pow(dy,2.0) + pow(dz,2.0));

        d_u = 2*(bcs.get_u(bc)- src.get_u(neighbour)) ;
        d_v = 2*(bcs.get_v(bc) - src.get_v(neighbour));
        d_w = 2*(bcs.get_w(bc) - src.get_w(neighbour));
        d_rho = 2*(bcs.get_rho(bc) - src.get_rho(neighbour)) ;

        n_i = mesh.get_face_i(i_face);
        n_j = mesh.get_face_j(i_face);
        n_k = mesh.get_face_k(i_face);

        if( rho_type){


            rho[i1].x = d_rho/d_mag*n_i;
            rho[i1].y = d_rho/d_mag*n_j;
            rho[i1].z = d_rho/d_mag*n_k;


        }else{
                u[i1].x =d_u/d_mag*n_i;
                v[i1].x = d_v/d_mag*n_i;
                w[i1].x = d_w/d_mag*n_i;

                u[i1].y = d_u/d_mag*n_j;
                v[i1].y = d_v/d_mag*n_j;
                w[i1].y = d_w/d_mag*n_j;

                u[i1].z = d_u/d_mag*n_k;
                v[i1].z = d_v/d_mag*n_k;
                w[i1].z = d_w/d_mag*n_k;
            }




    }

