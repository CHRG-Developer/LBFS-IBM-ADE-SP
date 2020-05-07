#include "tecplot_output.h"
#include "TECIO.h"
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include "post_processing.h"

#include "lagrangian_object.h"
#include <string>



tecplot_output::tecplot_output(){

}

tecplot_output::tecplot_output(global_variables &globals, Mesh &Mesh, Solution &Soln,
                             Boundary_Conditions &bcs, int fileType_e, double timestamp,
                             post_processing &pp)
{
    //ctor


    std::string output_location;
    std::string reynolds_text;
    reynolds_text = globals.reynolds_number;
    std::string zone_name;
    std::stringstream ss;
    std::stringstream tt;
    int *valueLocation = nullptr;
    INTEGER4 strandID;
    tt << timestamp;

    if( fileType_e == 1){
        output_location = globals.output_file + "/plt/grid.szplt";
        valueLocation = new int[3];
        for(int i = 0; i < 3;i++){
               valueLocation[i] = 1;
        }
        strandID  = 0;   /* StaticZone */

    }else{
        output_location = globals.output_file + "/plt/" + tt.str() +".szplt";
        valueLocation = new int[8];
        for(int i = 0; i < 8;i++){
               valueLocation[i] = 0;
        }
        strandID  = 1;
    }



   //enum fileType_e { FULL = 0, GRID = 1, SOLUTION = 2 };

     double *nodx, *nody, *z, *p,*u, *v,*y ,*x, *u_err, *u_exact ,*w, *vort, *st;
    int *connectivity;
    double solTime;
    INTEGER4 debug, i, j, k, dIsDouble, vIsDouble, zoneType,  parentZn, isBlock;
    INTEGER4 iCellMax, jCellMax, kCellMax, nFConns, fNMode, shrConn, fileType;


    INTEGER4 fileFormat; // 0 == PLT, 1 == SZPLT
    fileFormat = 1;

    INTEGER4 nNodes, nCells, nFaces, connectivityCount, index;
    int XDIM, YDIM, ZDIM; // nodes
    int t,r;
    int n_ghost;  // number of host cells;

    n_ghost = 1;

    if(globals.testcase == 3){
        n_ghost = 2;
    }

    XDIM = Mesh.get_num_x() +1 -2*n_ghost;
    YDIM = Mesh.get_num_y() +1 -2*n_ghost;
    ZDIM = 2;

    debug     = 1;
    vIsDouble = 0;
    dIsDouble = 0;
    nNodes = XDIM * YDIM * ZDIM;
    nCells = (XDIM - 1) * (YDIM - 1) * (ZDIM - 1);
    nFaces = 6; /* Not used */
    zoneType  = 5;      /* Brick */
    solTime   = timestamp;

    parentZn  = 0;      /* No Parent */
    isBlock   = 1;      /* Block */
    iCellMax  = 0;
    jCellMax  = 0;
    kCellMax  = 0;
    nFConns   = 0;
    fNMode    = 0;
    shrConn   = 0;
    fileType  = fileType_e;
    ss << nCells;
    zone_name = ss.str();
    /*
     * Open the file and write the tecplot datafile
     * header information
     */
     if( fileType_e == 1){
        i = TECINI142((char*) "Couette Flow" ,
                  (char*)"nodx nody z",
                  (char*) output_location.c_str(),
                  (char*) ".",
                  &fileFormat,
                  &fileType,
                  &debug,
                  &vIsDouble);


     }else{
     i = TECINI142((char*) "Couette Flow" ,
                  (char*)"p u v w x y vort stf",
                  (char*) output_location.c_str(),
                  (char*) ".",
                  &fileFormat,
                  &fileType,
                  &debug,
                  &vIsDouble);

     }

    i = TECAUXSTR142("Re" ,  reynolds_text.c_str());

     if( fileType_e == 1){
        nodx  = (double*)calloc(nNodes , sizeof(double));
        nody  = (double*)calloc(nNodes , sizeof(double));
        z  = (double*)calloc(nNodes , sizeof(double));
        t = 0;
        r=0;
        for (k = 0; k < ZDIM; k++){
            r=0;
            for (j = 0; j < Mesh.get_num_y() ; j++){
                for (i = 0; i < Mesh.get_num_x() ; i++)
                {
                    if( i != 0 && j != 0 && i != (n_ghost-1) && j != (n_ghost-1)
                    && i != (Mesh.get_num_x() - (n_ghost-1))
                     && j != (Mesh.get_num_y() - (n_ghost-1))){
                        nodx[t] = (double)(Mesh.get_west_x(r));
                        nody[t] = (double)(Mesh.get_south_y(r));
                        z[t] = (double)(k + 1);
                        t++;
                   }
                    r++;

                }
            }
        }

        connectivityCount = 8 * nCells;
        connectivity = (INTEGER4*)malloc(connectivityCount * sizeof(INTEGER4));
        for (k = 0; k < ZDIM - 1; k++)
            for (j = 0; j < YDIM - 1; j++)
                for (i = 0; i < XDIM - 1; i++)
                {
                    index = ((k * (YDIM - 1) + j) * (XDIM - 1) + i) * 8;
                    connectivity[index] = (k * YDIM + j) * XDIM + i +1;
                    connectivity[index + 1] = connectivity[index] + 1;
                    connectivity[index + 2] = connectivity[index] + XDIM + 1;
                    connectivity[index + 3] = connectivity[index] + XDIM;
                    connectivity[index + 4] = connectivity[index] + XDIM * YDIM;
                    connectivity[index + 5] = connectivity[index + 1] + XDIM * YDIM;
                    connectivity[index + 6] = connectivity[index + 2] + XDIM * YDIM;
                    connectivity[index + 7] = connectivity[index + 3] + XDIM * YDIM;
                }

     }
     else{
         p  = (double*)calloc(nCells , sizeof(double));
        u = (double*)calloc(nCells , sizeof(double));
        v = (double*)calloc(nCells , sizeof(double));
        w = (double*)calloc(nCells , sizeof(double));
        y = (double*)calloc(nCells , sizeof(double));
        x = (double*)calloc(nCells , sizeof(double));
        vort = (double*)calloc(nCells , sizeof(double));
        st = (double*)calloc(nCells , sizeof(double));
        t=0;
        for (i = 0; i < Mesh.get_total_cells(); ++i){
        if( bcs.get_bc(i) == false){
            p[t] = (double) Soln.get_rho(i);
            u[t] = (double) Soln.get_u(i)/globals.max_velocity;
            v[t] = (double) Soln.get_v(i)/globals.max_velocity;
            w[t] = (double) 0.0;
            x[t] = (double) Mesh.get_centroid_x(i)/Mesh.get_X();
            y[t] = (double) Mesh.get_centroid_y(i)/Mesh.get_Y();
            vort[t] = (double) pp.vorticity[i];
            st[t] = (double) pp.streamfunction[i];
            t++;
        }

    }

     }



    t = 0;



       /*
     * Write the zone header information.
     */
     std::cout << zone_name.c_str() << std::endl;
   i = TECZNE142(
                 //(char*) zone_name.c_str(),
                 zone_name.c_str(),
                  &zoneType,
                  &nNodes,
                  &nCells,
                  &nFaces,
                  &iCellMax,
                  &jCellMax,
                  &kCellMax,
                  &solTime,
                  &strandID,
                  &parentZn,
                  &isBlock,
                  &nFConns,
                  &fNMode,
                  0,              /* TotalNumFaceNodes */
                  0,              /* NumConnectedBoundaryFaces */
                  0,              /* TotalNumBoundaryConnections */
                  NULL,           /* PassiveVarList */
                  valueLocation,  /* ValueLocation = Nodal */
                  NULL,           /* SharVarFromZone */
                  &shrConn);
/*
     * Write out the field data.
     */

   if( fileType_e == 1){
        i = TECDAT142(&nNodes, nodx, &dIsDouble);
        i = TECDAT142(&nNodes, nody, &dIsDouble);
        i = TECDAT142(&nNodes, z, &dIsDouble);
       i = TECNODE142(&connectivityCount, connectivity);
        free(connectivity);
        free(nodx);
        free(nody);
        free(z);


   }else{
        i = TECDAT142(&nCells, p, &dIsDouble);
        i = TECDAT142(&nCells, u, &dIsDouble);
        i = TECDAT142(&nCells, v, &dIsDouble);
        i = TECDAT142(&nCells, w, &dIsDouble);
        i = TECDAT142(&nCells, x, &dIsDouble);
        i = TECDAT142(&nCells, y, &dIsDouble);
        i = TECDAT142(&nCells, vort, &dIsDouble);
        i = TECDAT142(&nCells, st, &dIsDouble);

        free(p);
        free(u);
        free(v);
        free(x);
        free(y);
        free(w);
        free(vort);
        free(st);

   }
    delete [] valueLocation;
    valueLocation= NULL;

    i = TECEND142();
}


void tecplot_output::tecplot_output_unstructured_soln(global_variables &globals, unstructured_mesh &Mesh, Solution &Soln,
                             Boundary_Conditions &bcs, int time,
                             post_processing &pp, Solution &residual, double *local_delta_t, double *visc_factor)
{
    //ctor

    int fileType_e = 2;
    std::string output_location;
    std::string reynolds_text;
    reynolds_text = globals.reynolds_number;
    std::string zone_name;
    std::stringstream ss;
    std::stringstream tt;
    int *valueLocation = nullptr;
    INTEGER4 strandID;
    tt << time;

    output_location = globals.output_file + "/plt/" + tt.str() +".plt";
    valueLocation = new int[16];
    for(int i = 0; i < 16;i++){
           valueLocation[i] = 0;
    }
    strandID  = 1;




   //enum fileType_e { FULL = 0, GRID = 1, SOLUTION = 2 };

     double *nodx, *nody, *z, *p,*u, *v,*y ,*x, *u_err, *u_exact ,*w, *vort, *st, *res_rho, *res_u,*res_v, *res_w, *dt, *visc;
    int *connectivity;
    double solTime;
    INTEGER4 debug, i, j, k, dIsDouble, vIsDouble, zoneType,  parentZn, isBlock;
    INTEGER4 iCellMax, jCellMax, kCellMax, nFConns, fNMode, shrConn, fileType;


    INTEGER4 fileFormat; // 0 == PLT, 1 == SZPLT
    fileFormat = 0;

    INTEGER4 nNodes, nCells, nFaces, connectivityCount, index;

    int t,r;


    debug     = 1;
    vIsDouble = 1;
    dIsDouble = 1;
    nNodes = Mesh.get_n_vertices();
    nCells = Mesh.get_n_cells();
    nFaces = Mesh.get_n_faces(); /* Not used */
    zoneType  = 7;      /* polyhedral */
    solTime   = time;

    parentZn  = 0;      /* No Parent */
    isBlock   = 1;      /* Block */
    iCellMax  = 0;
    jCellMax  = 0;
    kCellMax  = 0;
    nFConns   = 0;
    fNMode    = 0;
    shrConn   = 0;
    fileType  = fileType_e;
    ss << nCells;
    zone_name = ss.str();

    INTEGER4 NumFaceNodes = nFaces*4;
    INTEGER4 NumBConns  = 0;    /* No Boundary Connections */
    INTEGER4 NumBItems = 0;     /* No Boundary Items */

    /*
     * Open the file and write the tecplot datafile
     * header information
     */

     i = TECINI142((char*) "LBFS" ,
                  (char*)"p u v w x y z vort stf res_rho res_u res_v res_w dt visc",
                  (char*) output_location.c_str(),
                  (char*) ".",
                  &fileFormat,
                  &fileType,
                  &debug,
                  &vIsDouble);


    i = TECAUXSTR142("Re" ,  reynolds_text.c_str());

        p  = (double*)calloc(nCells , sizeof(double));
        u = (double*)calloc(nCells , sizeof(double));
        v = (double*)calloc(nCells , sizeof(double));
        w = (double*)calloc(nCells , sizeof(double));
        y = (double*)calloc(nCells , sizeof(double));
        x = (double*)calloc(nCells , sizeof(double));
        z = (double*)calloc(nCells , sizeof(double));
        vort = (double*)calloc(nCells , sizeof(double));
        st = (double*)calloc(nCells , sizeof(double));
		res_rho = (double*)calloc(nCells, sizeof(double));
		res_u = (double*)calloc(nCells, sizeof(double));
		res_v = (double*)calloc(nCells, sizeof(double));
		res_w = (double*)calloc(nCells, sizeof(double));
		dt = (double*)calloc(nCells, sizeof(double));
		visc = (double*)calloc(nCells, sizeof(double));

        t=0;
        for (i = 0; i < Mesh.get_n_cells(); ++i){

            p[t] = (double) Soln.get_rho(i);
            u[t] = (double) Soln.get_u(i) / globals.max_velocity;
            v[t] = (double) Soln.get_v(i) /globals.max_velocity;
            w[t] = (double) Soln.get_w(i) / globals.max_velocity;
            x[t] = (double) Mesh.get_centroid_x(i) /Mesh.get_X();
            y[t] = (double) Mesh.get_centroid_y(i) / Mesh.get_X();
            z[t] = (double) Mesh.get_centroid_z(i) / Mesh.get_X();
            vort[t] = (double) pp.vorticity[i];
            st[t] = (double) pp.streamfunction[i];
			res_rho[t] = (double) residual.get_rho(i);
			res_u[t] = (double) residual.get_u(i);
			res_v[t] = (double) residual.get_v(i);
			res_w[t] = (double) residual.get_w(i);
			dt[t] = (double)local_delta_t[i];
			visc[t] = (double)visc_factor[i];
			t++;
    }

    t = 0;

       /*
     * Write the zone header information.
     */
     std::cout << zone_name.c_str() << std::endl;
   i = TECZNE142(
                 //(char*) zone_name.c_str(),
                 zone_name.c_str(),
                  &zoneType,
                  &nNodes,
                  &nCells,
                  &nFaces,
                  &iCellMax,
                  &jCellMax,
                  &kCellMax,
                  &solTime,
                  &strandID,
                  &parentZn,
                  &isBlock,
                  &nFConns,
                  &fNMode,
                  &NumFaceNodes,              /* TotalNumFaceNodes */
                  &NumBConns,              /* NumConnectedBoundaryFaces */
                  &NumBItems,              /* TotalNumBoundaryConnections */
                  NULL,           /* PassiveVarList */
                  valueLocation,  /* ValueLocation = Nodal */
                  NULL,           /* SharVarFromZone */
                  &shrConn);
/*
     * Write out the field data.
     */

    i = TECDAT142(&nCells, p, &dIsDouble);
    i = TECDAT142(&nCells, u, &dIsDouble);
    i = TECDAT142(&nCells, v, &dIsDouble);
    i = TECDAT142(&nCells, w, &dIsDouble);
    i = TECDAT142(&nCells, x, &dIsDouble);
    i = TECDAT142(&nCells, y, &dIsDouble);
    i = TECDAT142(&nCells, z, &dIsDouble);
    i = TECDAT142(&nCells, vort, &dIsDouble);
    i = TECDAT142(&nCells, st, &dIsDouble);
	i = TECDAT142(&nCells, res_rho, &dIsDouble);
	i = TECDAT142(&nCells, res_u, &dIsDouble);
	i = TECDAT142(&nCells, res_v, &dIsDouble);
	i = TECDAT142(&nCells, res_w, &dIsDouble);
	i = TECDAT142(&nCells, dt, &dIsDouble);
	i = TECDAT142(&nCells, visc, &dIsDouble);

    free(p);
    free(u);
    free(v);
    free(x);
    free(y);
    free(z);
    free(w);
    free(vort);
    free(st);
	free(res_rho);
	free(res_u);
	free(res_v);
	free(res_w);
	free(dt);
	free(visc);


    delete [] valueLocation;
    valueLocation= NULL;

    i = TECEND142();
}





tecplot_output::~tecplot_output()
{
    //dtor
}




void tecplot_output::tecplot_output_lagrangian_object(lagrangian_object &object, global_variables &globals, domain_geometry &geometry, double timestamp)
{
	//ctor


	std::string output_location;
	std::string reynolds_text;
	reynolds_text = globals.reynolds_number;
	std::string zone_name;
	std::stringstream ss;
	std::stringstream tt;
	int *valueLocation = nullptr;
	INTEGER4 strandID;

	tt << timestamp;
	double X = geometry.X;
	double Y = geometry.Y;
	double Z = geometry.Z;
	double *nodx, *nody, *nodz, *vel_x, *vel_y, *vel_z, *force_x, *force_y, *force_z;
	int *connectivity;
	double solTime;
	INTEGER4 debug, i, j, k, dIsDouble, vIsDouble, zoneType, parentZn, isBlock;
	INTEGER4 iCellMax, jCellMax, kCellMax, nFConns, fNMode, shrConn, fileType;

	INTEGER4 nNodes, nCells, nFaces, connectivityCount, index;

	int t, r;
	int n_ghost;  // number of host cells;

	n_ghost = 1;

	if (globals.testcase == 3) {
		n_ghost = 2;
	}

	output_location = globals.output_file +  object.name + ".plt";
	valueLocation = new int[9];

	for (int i = 0; i < 9; i++) {
		valueLocation[i] = 1;
	}
	strandID = 5;   /* StaticZone */

   //enum fileType_e { FULL = 0, GRID = 1, SOLUTION = 2 };
	int fileType_e = 0;

	INTEGER4 fileFormat; // 0 == PLT, 1 == SZPLT
	fileFormat = 0;


	debug = 1;
	vIsDouble = 1;
	dIsDouble = 1;
	//num nodes and num_cells will be changed in the future
	nNodes = object.num_nodes ;
	nCells = object.radial_nodes*(object.depth_nodes - 1);
	nFaces = object.radial_nodes * object.depth_nodes + object.radial_nodes*(object.depth_nodes-1) ; /* Not used */
	zoneType = 6;      /* FEPolygon */
	solTime = timestamp;

	parentZn = 0;      /* No Parent */
	isBlock = 1;      /* Block */
	iCellMax = 0;
	jCellMax = 0;
	kCellMax = 0;
	nFConns = 0;
	fNMode = 1;
	shrConn = 0;
	fileType = fileType_e;
	ss << nCells;
	zone_name = ss.str();

	/* The number of face nodes in the zone. This example creates
	* a zone with a single pyramidal cell. This cell has four
	* triangular faces and one rectangular face, yielding a total
	* of 16 face nodes.
	*/
	INTEGER4 NumFaceNodes = nFaces * 2;
	INTEGER4 NumBConns = 0;    /* No Boundary Connections */
	INTEGER4 NumBItems = 0;     /* No Boundary Items */


	/*
	 * Open the file and write the tecplot datafile
	 * header information
	 */
	i = TECINI142((char*)object.name.c_str(),
		(char*)"nodx nody nodz vel_x vel_y vel_z force_x force_y force_z",
		(char*)output_location.c_str(),
		(char*) ".",
		&fileFormat,
		&fileType,
		&debug,
		&vIsDouble);




	i = TECAUXSTR142("Re", reynolds_text.c_str());


	nodx = (double*)calloc(nNodes, sizeof(double));
	nody = (double*)calloc(nNodes, sizeof(double));
	nodz = (double*)calloc(nNodes, sizeof(double));
	vel_x = (double*)calloc(nNodes, sizeof(double));
	vel_y = (double*)calloc(nNodes, sizeof(double));
	vel_z = (double*)calloc(nNodes, sizeof(double));
	force_x = (double*)calloc(nNodes, sizeof(double));
	force_y = (double*)calloc(nNodes, sizeof(double));
	force_z = (double*)calloc(nNodes, sizeof(double));

	t = 0;
	r = 0;

	std::string filename;
	std::ofstream globals_txt, globals_txt1;
	std::string globals_file, globals_file1;
	output_location = globals.output_file;
	filename = globals.simulation_name;
	globals_file = output_location + "/trial.txt";
	globals_file1 = output_location + "/trial1.txt";
	globals_txt.open(globals_file.c_str(), ios::out);


	for (k = 0; k < object.num_nodes; k++) {

		nodx[t] = (double)(object.node_x[k]) / X;
		nody[t] = (double)(object.node_y[k]) / Y;
		nodz[t] = (double)(object.node_z[k]) / Z; // temp add z coordinates
		vel_x[t] = (double)(object.node_vel_x[k]) / globals.max_velocity;
		vel_y[t] = (double)(object.node_vel_y[k]) / globals.max_velocity;
		vel_z[t] = (double)(object.node_vel_z[k]) / globals.max_velocity;

		force_x[t] = (double)(object.node_force_x[k]) / globals.max_velocity;
		force_y[t] = (double)(object.node_force_y[k]) / globals.max_velocity;
		force_z[t] = (double)(object.node_force_z[k]) / globals.max_velocity;

		globals_txt << t << "," << nodx[t] << "," << nody[t] << "," << nodz[t] << std::endl;
		t++;

	}

	globals_txt.close();

	/*
* Write the zone header information.
*/
	std::cout << zone_name.c_str() << std::endl;
	i = TECZNE142(
		//(char*) zone_name.c_str(),
		zone_name.c_str(),
		&zoneType,
		&nNodes,
		&nCells,
		&nFaces,
		&iCellMax,
		&jCellMax,
		&kCellMax,
		&solTime,
		&strandID,
		&parentZn,
		&isBlock,
		&nFConns,
		&fNMode,
		&NumFaceNodes,              /* TotalNumFaceNodes */
		&NumBConns,              /* NumConnectedBoundaryFaces */
		&NumBItems,              /* TotalNumBoundaryConnections */
		NULL,           /* PassiveVarList */
		valueLocation,  /* ValueLocation = Nodal */
		NULL,           /* SharVarFromZone */
		&shrConn);
	/*
		 * Write out the field data.
		 */

	i = TECDAT142(&nNodes, nodx, &dIsDouble);
	i = TECDAT142(&nNodes, nody, &dIsDouble);
	i = TECDAT142(&nNodes, nodz, &dIsDouble);

	i = TECDAT142(&nNodes, vel_x, &dIsDouble);
	i = TECDAT142(&nNodes, vel_y, &dIsDouble);
	i = TECDAT142(&nNodes, vel_z, &dIsDouble);

	i = TECDAT142(&nNodes, force_x, &dIsDouble);
	i = TECDAT142(&nNodes, force_y, &dIsDouble);
	i = TECDAT142(&nNodes, force_z, &dIsDouble);


	INTEGER4 *FaceNodes = new INTEGER4[NumFaceNodes];
	int q = 0;


	//radial faces
	for (int j = 0; j < object.depth_nodes; j++) {
		for (int k = 1; k < object.radial_nodes; k++) {

			FaceNodes[q] = object.radial_nodes*j + k;
			q++;
			FaceNodes[q] = object.radial_nodes*j + k + 1;
			q++;
		}
		FaceNodes[q] = object.radial_nodes * (j + 1);
		q++;
		FaceNodes[q] = 1 + object.radial_nodes*j;
		q++;
	}

	//axial faces

	for (int j = 0; j < object.depth_nodes-1; j++) {
		for (int k = 0; k < object.radial_nodes; k++) {

			FaceNodes[q] = object.radial_nodes * j + k+ 1;
			q++;
			FaceNodes[q] = object.radial_nodes * (j + 1) +  k + 1;
			q++;
		}

	}


	INTEGER4 *FaceLeftElems = new INTEGER4[nFaces];
	INTEGER4 *FaceRightElems = new INTEGER4[nFaces];
	int NF = 0;

	//first border of radial faces
	for (INTEGER4 ii = 0; ii < object.radial_nodes; ii++)
	{
		FaceLeftElems[NF] = 0;

		FaceRightElems[NF] = ii +1;
		NF++;
	}

	/// interior
	for (int j = 1; j < object.depth_nodes-1; j++) {
		for (INTEGER4 ii = 0; ii < object.radial_nodes; ii++)
		{
			FaceLeftElems[NF] = object.radial_nodes *(j-1) + ii + 1;

			FaceRightElems[NF] = object.radial_nodes *(j ) + ii + 1;
			NF++;
		}
	}

	//right border
	for (INTEGER4 ii = 0; ii < object.radial_nodes; ii++)
	{
		FaceLeftElems[NF] = object.radial_nodes *(object.depth_nodes - 2) + ii + 1;

		FaceRightElems[NF] = 0;
		NF++;
	}

	//axial elements

	for (int j = 0; j < object.depth_nodes-1; j++) {

		//very last cell
		FaceLeftElems[NF] = object.radial_nodes *(j)+1;
		FaceRightElems[NF] = object.radial_nodes *(j + 1);
		NF++;

		for (INTEGER4 ii = 1; ii < object.radial_nodes; ii++)
		{
			FaceLeftElems[NF] = object.radial_nodes *(j)+ii +1 ;

			FaceRightElems[NF] = object.radial_nodes *(j)+ii ;
			NF++;
		}


	}

	i = tecpolyface142(&nFaces, NULL,  FaceNodes, FaceLeftElems, FaceRightElems);

	i = TECEND142();

	delete FaceNodes;
	delete FaceLeftElems;
	delete FaceRightElems;

	free(nodx);
	free(nody);
	free(nodz);
	free(vel_x);
	free(vel_y);
	free(vel_z);
	free(force_x);
	free(force_y);
	free(force_z);



	delete[] valueLocation;
	valueLocation = NULL;



}


void tecplot_output::spring_network_output(std::string fileName, bool Ctn, lagrangian_object &object) {

	//coordinates output.
	std::cout << "Preparing for data output... ";

	std::ofstream outputfile;
	std:string output_location;

	output_location = fileName + "/plt/object/" + object.name + ".dat";

	if (Ctn)
		outputfile.open(output_location, std::ofstream::out | std::ofstream::app);
	else
		outputfile.open(output_location);

	outputfile << "VARIABLES = \"x\" \"y\" \"z\" "
		<< "\"i_force\" \"j_force\" \"k_force\" \"curvature\" \"index\""
		<< std::endl;
	outputfile << "ZONE Nodes= " << (object.num_nodes) << " Elements= " << (object.num_tets) << " "
		<< "F=FEPOINT ET=TRIANGLE \n"
		<< "DT=(DOUBLE DOUBLE DOUBLE "
		<< "DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE "
		<< ")" << std::endl;

	for (int j = 0; j < object.num_nodes; j++) {


		outputfile << object.node_x[j] << "\t";
		outputfile << object.node_y[j] << "\t";
		outputfile << object.node_z[j] << "\t";

		outputfile << object.node_force_x[j] << "\t";
		outputfile << object.node_force_y[j] << "\t";
		outputfile << object.node_force_z[j] << "\t";

		outputfile << object.node_vel_x[j] << "\t";
		outputfile << j << "\t";

		outputfile << std::endl;
	}

	for (int i = 0; i < object.num_tets; i++) {
		for (int j = 0; j < 3; j++) {
			outputfile << object.tet_connectivity[3 * i + j] + 1 << "\t";

		}
		outputfile << std::endl;
	}

	outputfile.close();

	std::cout << "Data output successful! " << std::endl;

	return;
}




void tecplot_output::spring_network_output_gpu(std::string fileName, bool Ctn, double * obj_vel_x, double * obj_vel_y, double * obj_vel_z, double * obj_x, double * obj_y, double * obj_z,
	double * obj_force_x, double * obj_force_y, double * obj_force_z, std::string name,  int time, int num_nodes, int num_tets, int * tet_connectivity, double * nodal_area, bool final, double * curvature) {

	//coordinates output.
	std::cout << "Preparing for data output... ";

	std::ofstream outputfile;
	std:string output_location;
	std::stringstream tt;

	tt << time;
	if (final) {
		output_location = fileName + "/" + name  + ".dat";
	}
	else {
		output_location = fileName + "/plt/object/" + name + tt.str() + ".dat";
	}


	if (Ctn)
		outputfile.open(output_location, std::ofstream::out | std::ofstream::app);
	else
		outputfile.open(output_location);

	outputfile << "VARIABLES = \"nodx\" \"nody\" \"nodz\" "
		<< "\"i_force\" \"j_force\" \"k_force\" \"vel_x\" \"vel_y\" \"vel_z\" \"nodal_area\" \"index\" \"curvature\""
		<< std::endl;
	outputfile << "ZONE Nodes= " << (num_nodes) << " Elements= " << (num_tets) << " "
		<< "F=FEPOINT ET=TRIANGLE \n"

		<< "SOLUTIONTIME = " << tt.str() << " \n"
		<< "DT=(DOUBLE DOUBLE DOUBLE "
		<< "DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE"
		<< ")" << std::endl;

	for (int j = 0; j < num_nodes; j++) {


		outputfile << obj_x[j] << "\t";
		outputfile << obj_y[j] << "\t";
		outputfile << obj_z[j] << "\t";

		outputfile << obj_force_x[j] << "\t";
		outputfile << obj_force_y[j] << "\t";
		outputfile << obj_force_z[j] << "\t";

		outputfile << obj_vel_x[j] << "\t";
		outputfile << obj_vel_y[j] << "\t";
		outputfile << obj_vel_z[j] << "\t";

		outputfile << nodal_area[j] << "\t";
		outputfile << j << "\t";
		outputfile << curvature[j] << "\t";
		outputfile << std::endl;
	}

	for (int i = 0; i < num_tets ; i++) {
		for (int j = 0; j < 3; j++) {
			outputfile << tet_connectivity[3 * i + j] + 1 << "\t";
			outputfile << std::endl;
		}

	}

	outputfile.close();

	std::cout << "Data output successful! " << std::endl;

	return;
}

//
//
//void tecplot_output::spring_network_output_gpu_post_processing(std::string fileName, bool Ctn, double * obj_vel_x, double * obj_vel_y, double * obj_vel_z, double * obj_x, double * obj_y, double * obj_z,
//	double * obj_force_x, double * obj_force_y, double * obj_force_z, std::string name, int time, int num_nodes, int num_tets, int * tet_connectivity, double * nodal_area, bool final) {
//
//	//coordinates output.
//	std::cout << "Preparing for data output... ";
//
//	std::ofstream outputfile;
//std:string output_location;
//	std::stringstream tt;
//
//	tt << time;
//	if (final) {
//		output_location = fileName + "/" + name + ".dat";
//	}
//	else {
//		output_location = fileName + "/plt/object/" + name + tt.str() + ".dat";
//	}
//
//
//	if (Ctn)
//		outputfile.open(output_location, std::ofstream::out | std::ofstream::app);
//	else
//		outputfile.open(output_location);
//
//	outputfile << "VARIABLES = \"nodx\" \"nody\" \"nodz\" "
//		<< "\"i_force\" \"j_force\" \"k_force\" \"vel_x\" \"vel_y\" \"vel_z\" \"nodal_area\" "
//		<< "\"curvature\" "
//		<< "\"dilate\" \"shear\" "
//		<< std::endl;
//	outputfile << "ZONE Nodes= " << (num_nodes) << " Elements= " << (num_tets) << " "
//		<< "ZONETYPE=FETRIANGLE \n"
//		<< "DATAPACKING=BLOCK \n"
//		<< "VARLOCATION=([8-9]=CELLCENTERED) \n"
//
//		<< "SOLUTIONTIME = " << tt.str() << " \n"
//		<< "DT=(DOUBLE DOUBLE DOUBLE "
//		<< "DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE"
//		<< ")" << std::endl;
//
//	for (int j = 0; j < num_nodes; j++) {
//
//
//		outputfile << obj_x[j] << "\t";
//		outputfile << obj_y[j] << "\t";
//		outputfile << obj_z[j] << "\t";
//
//		outputfile << obj_force_x[j] << "\t";
//		outputfile << obj_force_y[j] << "\t";
//		outputfile << obj_force_z[j] << "\t";
//
//		outputfile << obj_vel_x[j] << "\t";
//		outputfile << obj_vel_y[j] << "\t";
//		outputfile << obj_vel_z[j] << "\t";
//
//		outputfile << nodal_area[j] << "\t";
//
//		outputfile << std::endl;
//	}
//
//	for (int i = 0; i < num_tets; i++) {
//		for (int j = 0; j < 3; j++) {
//			outputfile << tet_connectivity[3 * i + j] + 1 << "\t";
//			outputfile << std::endl;
//		}
//
//	}
//
//	outputfile.close();
//
//	std::cout << "Data output successful! " << std::endl;
//
//	return;
//}

void tecplot_output::tecplot_output_lagrangian_object_gpu(double * obj_vel_x, double * obj_vel_y, double * obj_vel_z, double * obj_x, double * obj_y, double * obj_z,
	double * obj_force_x, double * obj_force_y, double * obj_force_z,
	global_variables &globals, domain_geometry &geometry, double timestamp, std::string obj_name, int obj_nodes, int depth_nodes, int radial_nodes)
{
	//ctor


	std::string output_location;
	std::string reynolds_text;
	reynolds_text = globals.reynolds_number;
	std::string zone_name;
	std::stringstream ss;
	std::stringstream tt;
	int *valueLocation = nullptr;
	INTEGER4 strandID;

	tt << timestamp;
	double X = geometry.X;
	double Y = geometry.Y;
	double Z = geometry.Z;
	double *nodx, *nody, *nodz, *vel_x, *vel_y, *vel_z, *force_x, *force_y, *force_z;
	int *connectivity;
	double solTime;
	INTEGER4 debug, i, j, k, dIsDouble, vIsDouble, zoneType, parentZn, isBlock;
	INTEGER4 iCellMax, jCellMax, kCellMax, nFConns, fNMode, shrConn, fileType;

	INTEGER4 nNodes, nCells, nFaces, connectivityCount, index;

	int t, r;
	int n_ghost;  // number of host cells;

	n_ghost = 1;

	if (globals.testcase == 3) {
		n_ghost = 2;
	}

	output_location = globals.output_file + "/plt/object/" + obj_name + tt.str() + ".plt";
	valueLocation = new int[9];

	for (int i = 0; i < 9; i++) {
		valueLocation[i] = 1;
	}
	strandID = 5;   /* StaticZone */

   //enum fileType_e { FULL = 0, GRID = 1, SOLUTION = 2 };
	int fileType_e = 0;

	INTEGER4 fileFormat; // 0 == PLT, 1 == SZPLT
	fileFormat = 0;


	debug = 1;
	vIsDouble = 1;
	dIsDouble = 1;
	//num nodes and num_cells will be changed in the future
	nNodes = obj_nodes;
	nCells = radial_nodes*(depth_nodes - 1);
	nFaces = radial_nodes * depth_nodes + radial_nodes*(depth_nodes - 1); /* Not used */
	zoneType = 6;      /* FEPolygon */
	solTime = timestamp;

	parentZn = 0;      /* No Parent */
	isBlock = 1;      /* Block */
	iCellMax = 0;
	jCellMax = 0;
	kCellMax = 0;
	nFConns = 0;
	fNMode = 1;
	shrConn = 0;
	fileType = fileType_e;
	ss << nCells;
	zone_name = ss.str();

	/* The number of face nodes in the zone. This example creates
	* a zone with a single pyramidal cell. This cell has four
	* triangular faces and one rectangular face, yielding a total
	* of 16 face nodes.
	*/
	INTEGER4 NumFaceNodes = nFaces * 2;
	INTEGER4 NumBConns = 0;    /* No Boundary Connections */
	INTEGER4 NumBItems = 0;     /* No Boundary Items */


	/*
	 * Open the file and write the tecplot datafile
	 * header information
	 */
	i = TECINI142((char*)obj_name.c_str(),
		(char*)"nodx nody nodz vel_x vel_y vel_z force_x force_y force_z",
		(char*)output_location.c_str(),
		(char*) ".",
		&fileFormat,
		&fileType,
		&debug,
		&vIsDouble);




	i = TECAUXSTR142("Re", reynolds_text.c_str());


	nodx = (double*)calloc(nNodes, sizeof(double));
	nody = (double*)calloc(nNodes, sizeof(double));
	nodz = (double*)calloc(nNodes, sizeof(double));
	vel_x = (double*)calloc(nNodes, sizeof(double));
	vel_y = (double*)calloc(nNodes, sizeof(double));
	vel_z = (double*)calloc(nNodes, sizeof(double));
	force_x = (double*)calloc(nNodes, sizeof(double));
	force_y = (double*)calloc(nNodes, sizeof(double));
	force_z = (double*)calloc(nNodes, sizeof(double));

	t = 0;
	r = 0;

	std::string filename;
	std::ofstream globals_txt, globals_txt1;
	std::string globals_file, globals_file1;
	output_location = globals.output_file;
	filename = globals.simulation_name;
	globals_file = output_location + "/trial.txt";
	globals_file1 = output_location + "/trial1.txt";
	globals_txt.open(globals_file.c_str(), ios::out);


	for (k = 0; k < obj_nodes; k++) {

		nodx[t] = (double)(obj_x[k]) / X;
		nody[t] = (double)(obj_y[k]) / Y;
		nodz[t] = (double)(obj_z[k]) / Z; // temp add z coordinates
		vel_x[t] = (double)(obj_vel_x[k]) ;
		vel_y[t] = (double)(obj_vel_y[k]) ;
		vel_z[t] = (double)(obj_vel_z[k]);

		force_x[t] = (double)(obj_force_x[k]) ;
		force_y[t] = (double)(obj_force_y[k]) ;
		force_z[t] = (double)(obj_force_z[k]) ;

		globals_txt << t << "," << nodx[t] << "," << nody[t] << "," << nodz[t] << std::endl;
		t++;

	}

	globals_txt.close();

	/*
* Write the zone header information.
*/
	std::cout << zone_name.c_str() << std::endl;
	i = TECZNE142(
		//(char*) zone_name.c_str(),
		zone_name.c_str(),
		&zoneType,
		&nNodes,
		&nCells,
		&nFaces,
		&iCellMax,
		&jCellMax,
		&kCellMax,
		&solTime,
		&strandID,
		&parentZn,
		&isBlock,
		&nFConns,
		&fNMode,
		&NumFaceNodes,              /* TotalNumFaceNodes */
		&NumBConns,              /* NumConnectedBoundaryFaces */
		&NumBItems,              /* TotalNumBoundaryConnections */
		NULL,           /* PassiveVarList */
		valueLocation,  /* ValueLocation = Nodal */
		NULL,           /* SharVarFromZone */
		&shrConn);
	/*
		 * Write out the field data.
		 */

	i = TECDAT142(&nNodes, nodx, &dIsDouble);
	i = TECDAT142(&nNodes, nody, &dIsDouble);
	i = TECDAT142(&nNodes, nodz, &dIsDouble);

	i = TECDAT142(&nNodes, vel_x, &dIsDouble);
	i = TECDAT142(&nNodes, vel_y, &dIsDouble);
	i = TECDAT142(&nNodes, vel_z, &dIsDouble);

	i = TECDAT142(&nNodes, force_x, &dIsDouble);
	i = TECDAT142(&nNodes, force_y, &dIsDouble);
	i = TECDAT142(&nNodes, force_z, &dIsDouble);


	INTEGER4 *FaceNodes = new INTEGER4[NumFaceNodes];
	int q = 0;


	//radial faces
	for (int j = 0; j < depth_nodes; j++) {
		for (int k = 1; k < radial_nodes; k++) {

			FaceNodes[q] = radial_nodes*j + k;
			q++;
			FaceNodes[q] = radial_nodes*j + k + 1;
			q++;
		}
		FaceNodes[q] = radial_nodes * (j + 1);
		q++;
		FaceNodes[q] = 1 + radial_nodes*j;
		q++;
	}

	//axial faces

	for (int j = 0; j < depth_nodes - 1; j++) {
		for (int k = 0; k < radial_nodes; k++) {

			FaceNodes[q] = radial_nodes * j + k + 1;
			q++;
			FaceNodes[q] = radial_nodes * (j + 1) + k + 1;
			q++;
		}

	}


	INTEGER4 *FaceLeftElems = new INTEGER4[nFaces];
	INTEGER4 *FaceRightElems = new INTEGER4[nFaces];
	int NF = 0;

	//first border of radial faces
	for (INTEGER4 ii = 0; ii < radial_nodes; ii++)
	{
		FaceLeftElems[NF] = 0;

		FaceRightElems[NF] = ii + 1;
		NF++;
	}

	/// interior
	for (int j = 1; j < depth_nodes - 1; j++) {
		for (INTEGER4 ii = 0; ii < radial_nodes; ii++)
		{
			FaceLeftElems[NF] = radial_nodes *(j - 1) + ii + 1;

			FaceRightElems[NF] = radial_nodes *(j)+ii + 1;
			NF++;
		}
	}

	//right border
	for (INTEGER4 ii = 0; ii < radial_nodes; ii++)
	{
		FaceLeftElems[NF] = radial_nodes *(depth_nodes - 2) + ii + 1;

		FaceRightElems[NF] = 0;
		NF++;
	}

	//axial elements

	for (int j = 0; j < depth_nodes - 1; j++) {

		//very last cell
		FaceLeftElems[NF] = radial_nodes *(j)+1;
		FaceRightElems[NF] = radial_nodes *(j + 1);
		NF++;

		for (INTEGER4 ii = 1; ii < radial_nodes; ii++)
		{
			FaceLeftElems[NF] = radial_nodes *(j)+ii + 1;

			FaceRightElems[NF] = radial_nodes *(j)+ii;
			NF++;
		}


	}

	i = tecpolyface142(&nFaces, NULL, FaceNodes, FaceLeftElems, FaceRightElems);

	i = TECEND142();

	delete FaceNodes;
	delete FaceLeftElems;
	delete FaceRightElems;

	free(nodx);
	free(nody);
	free(nodz);
	free(vel_x);
	free(vel_y);
	free(vel_z);
	free(force_x);
	free(force_y);
	free(force_z);



	delete[] valueLocation;
	valueLocation = NULL;



}
