#include "unstructured_mesh.h"

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
#include "vector_var.h"
#include <numeric>
#include <cmath>
#include <tuple>

using namespace std;


unstructured_mesh::unstructured_mesh(domain_geometry &domain, global_variables &globals) : Mesh( domain, globals)
{
    //ctor

    if( globals.mesh_type < 2){
        return;
    }

    int i;
    X = domain.X;
    Y = domain.Y;
    Z = domain.Z;
    dx = domain.dx; //dimensional form
    dy = domain.dy;
    dz = domain.dz;

    //declare variables

    // get_CGNS_mesh
    if(globals.import_format == "openfoam"){
		if (globals.mesh_type == 3) {
			import_openfoam_mesh(globals);
		}else {
			create_uniform_openfoam_mesh(globals, domain);
		}
		openfoam_connectivity(globals);
        tecplot_output_polyhedral_unstructured_grid(globals);
        i = i +1;



    }else if(globals.import_format == "CGNS"){
//        open_cgns_file(globals);
//
//        get_num_bc_elements(globals);
//        import_cgns_mesh();
//        import_cgns_bcs();
//        //Output Tecplot Grid
//        tecplot_output_unstructured_grid(globals);

    }


    i = i +1;


    // Transfer to internal data mapping

    initialise_mesh_variables();

    generate_internal_mesh(globals);

    if(globals.import_format == "openfoam"){
        import_openfoam_BC_tags(globals);
		std::cout << "get immersed boundary look ups" << endl;
		generate_immersed_boundary_look_up(globals,domain);
		std::cout << "complete" << endl;
    }else{
//        read_bc_tags(globals);
    }

    //perform calcs on face neighbour/centorid etc.
    output_mesh_to_text(globals);


    destruct_mesh_variables(globals);

    if(globals.import_format == "openfoam"){

    }else{
//          /* close CGNS file */
//        cg_close(index_file);
    }


}

unstructured_mesh::~unstructured_mesh()
{
    //dtor
    
    delete [](calcs_per_cell);
    calcs_per_cell = NULL;

    delete [](bc_neighbour);
    bc_neighbour = NULL;

     delete [](bc_face_x);
    bc_face_x = NULL;
    delete [](bc_face_y);
    bc_face_y = NULL;
     delete [](bc_face_z);
    bc_face_z = NULL;

    delete [](bc_face_direction);
    bc_face_direction = NULL;

    delete[] (total_cell_area);
    total_cell_area = NULL;

    delete[] (face_x);
    face_x = NULL;

    delete[] (face_y);
    face_y = NULL;

    delete[] (face_z);
    face_z = NULL;

    delete[] (face_area);
    face_area = NULL;

    delete[] (delta_t_face);
    delta_t_face = NULL;

    delete[] (owner);
    owner = NULL;

    delete[] (neighbour);
    neighbour = NULL;

   

}
void unstructured_mesh::openfoam_connectivity(global_variables &globals){

//    //sort nodes by xyz
//    std::sort(openfoam_map.begin(), openfoam_map.end() , doCompare(*this));

    int i,n;
    i = 0;
    int pos,count_zero;

    std::vector<double> x_vol,y_vol,z_vol;
    std::vector<int> vol_index(8);
    x_vol.push_back(0.0);
    y_vol.push_back(0.0);
    z_vol.push_back(0.0);

    double x_min, x_max, y_min,y_max, z_min, z_max;

    //volume connectivity

    for(int i =0 ; i < n_faces; i++){
        volume_connectivity.push_back(std::vector <int> ());
    }

    for(int i =0 ; i < n_faces; i++){

        volume_connectivity[owner[i]].push_back(i);

        if( i < n_neighbours){
            volume_connectivity[neighbour[i]].push_back(i);

        }
    }


}
void unstructured_mesh::generate_immersed_boundary_look_up(global_variables &globals, domain_geometry &domain) {

	//get boundary condition face
	// example inlet

	
	
	int t =0;
	int j;

	double temp_x, temp_y;

	for (int i = 0;  i < num_bc_cells; i++) {
		if (bc_types[i].compare(domain.IBM_face)  == 0 ) {

			j = owner[i + n_neighbours];

			//IBM_lookup.push_back(std::make_tuple(face_x[n_neighbours + i], face_y[n_neighbours + i], t));
			temp_x = ceil(centroid_x[j] / 0.25)*0.25;
			temp_y = ceil(centroid_y[j] / 0.25)*0.25;

			IBM_lookup.push_back(std::make_tuple(temp_x, temp_y, j));
			t++;
		}
	}

	
		
	//now sort into x-y rank to make binary sort possible for efficiencies later on
	sort(IBM_lookup.begin(), IBM_lookup.end(), [](std::tuple < double,double, int> const &t1, std::tuple < double, double, int> const &t2) 
		{
			return (std::get<0>(t1) > std::get<0>(t2) ||
				(fabs(std::get<0>(t1) -std::get<0>(t2)) < 0.0001 ) && (std::get<1>(t1) > std::get<1>(t2)) ||
				(fabs(std::get<0>(t1) - std::get<0>(t2)) < 0.0001) && (fabs(std::get<1>(t1) - std::get<1>(t2)) < 0.0001)

				) ; // or use a custom compare function
			// && (std::get<1>(t1) > std::get<1>(t2)
		}
	);

	for (int i = 0; i < IBM_lookup.size(); i++) {
		cout << get<0>(IBM_lookup[i]) << " "
			<< get<1>(IBM_lookup[i]) << " "
			<< get<2>(IBM_lookup[i]) << "\n";
	}


	//declare and pass values to plateau functions x, y, z

	//will need slight modification for force spreading etc. 
	
}

void unstructured_mesh::create_uniform_openfoam_mesh(global_variables &globals, domain_geometry &domain) {
	std::string input_file;
	std::ifstream points;
	std::string line;
	int t, i,j,k;
	int a, b, c, d;
	n_vertices = 0;

	//read in owners file

	//calculate number of cells

	int n_cells_x, n_cells_y, n_cells_z;

	n_cells_x = ceil(domain.X / domain.dt*0.5);
	n_cells_y = ceil(domain.Y / domain.dt*0.5);
	n_cells_z = ceil(domain.Z / domain.dt*0.5);

	double dx, dy, dz;

	dx = domain.X / n_cells_x;
	dy = domain.Y / n_cells_y;
	dz = domain.Z / n_cells_z;

	delta_h = min(dx, min(dy, dz)) * 0.5;
	delta_h = min(delta_h, 0.1);

	double sin_arg;
	n_vertices = (n_cells_x + 1) *(n_cells_y + 1) * (n_cells_z + 1);
	int n_vertices_x, n_vertices_y, n_vertices_z;

	n_vertices_x = n_cells_x + 1;
	n_vertices_y = n_cells_y + 1;
	n_vertices_z = n_cells_z + 1;


	x = new float[n_vertices];
	if (x == NULL) exit(1);
	y = new float[n_vertices];
	if (y == NULL) exit(1);
	z = new float[n_vertices];
	if (z == NULL) exit(1);

	
	plateau_x = new double[n_cells_x + 1];
	plateau_y = new double[n_cells_y + 1];
	plateau_z = new double[n_cells_z + 1];

	double test_i, test_f;

	t = 0;

	double fx,fy,fz;

	double ap, bp;  /// tuning parameters of the plateau function
	ap = domain.plateau_width;
	bp = domain.plateau_slope;
	double  plateau_sum_x, plateau_sum_y, plateau_sum_z;
	plateau_sum_x = 0.0;

	// need to get sum of x,y,z plateau lengths for normalising to actual lengths
	// could save run time by saving to array but not arsed with memory management
	for (int i = 1; i < n_cells_x + 1; i++) {
		fx = i - n_cells_x / 2.0 -0.5; 
		plateau_x[i] = (exp(bp*(fx - ap)) + 1) * (exp(bp*(-fx - ap)) + 1);
		plateau_sum_x = plateau_sum_x + plateau_x[i];
	}

	for (int i = 1; i < n_cells_y + 1; i++) {
		fx = i - n_cells_y / 2.0 -0.5;
		plateau_y[i] = (exp(bp*(fx - ap)) + 1) * (exp(bp*(-fx - ap)) + 1);
		plateau_sum_y = plateau_sum_y + plateau_y[i];
	}

	for (int i = 1; i < n_cells_z + 1; i++) {
		fx = i - n_cells_z / 2.0 -0.5;
		plateau_z[i] = (exp(bp*(fx - ap)) + 1) * (exp(bp*(-fx - ap)) + 1);
		plateau_sum_z = plateau_sum_z + plateau_z[i];
	}

	for (int i = 0; i < n_cells_x + 1; i++) {
		if (i == 0) {
			plateau_x[i] = domain.origin_x;
		}
		else {
			plateau_x[i] = plateau_x[i - 1] + plateau_x[i] * domain.X / plateau_sum_x;
		}
	}

	for (int i = 0; i < n_cells_y + 1; i++) {
		if (i == 0) {
			plateau_y[i] = domain.origin_y;
		}
		else {
			plateau_y[i] = plateau_y[i - 1] + plateau_y[i] * domain.Y / plateau_sum_y;
		}
	}

	for (int i = 0; i < n_cells_z + 1; i++) {
		if (i == 0) {
			plateau_z[i] = domain.origin_z;
		}
		else {
			plateau_z[i] = plateau_z[i - 1] + plateau_z[i] * domain.Z / plateau_sum_z;
		}
	}

	/*int test;
	test = findClosest(plateau_x, n_cells_x + 1, -5.0);
*/
	   
	// note n_node = n_cell +1 
	for (int k = 0; k < n_cells_z + 1; k++) {
		for (int j = 0; j < n_cells_y +1; j++) {
			for (int i = 0; i < n_cells_x +1; i++) {
				/**/
				
				if (globals.mesh_type == 6) {
					x[t] = domain.origin_x + i * dx;
					y[t] = domain.origin_y + j * dy;
					z[t] = domain.origin_z + k * dz;

				}
				else if (globals.mesh_type == 7) {
					x[t] = domain.origin_x + 0.5 *(1 - cos(i / (double)(n_cells_x)* globals.PI)) * domain.X;
					y[t] = domain.origin_y + 0.5 *(1 - cos(j / (double)(n_cells_y)* globals.PI)) * domain.Y;
					z[t] = domain.origin_z + k * dz;

				}
				else if (globals.mesh_type == 4) {
					x[t] = plateau_x[i];
					y[t] = plateau_y[j];
					z[t] = plateau_z[k];
				}
				else if (globals.mesh_type == 5) {
					x[t] = domain.origin_x + 0.5 *( 1- cos( i / (double) (n_cells_x)  * globals.PI ) ) * domain.X;
					y[t] = domain.origin_y + 0.5 *(1 - cos(j/ (double)(n_cells_y )  * globals.PI)) * domain.Y;
					z[t] = domain.origin_z + 0.5 *(1 - cos(k / (double)(n_cells_z )  * globals.PI)) * domain.Z;
				}
			t++;
			}
		}
	}

	t = t;

	n_cells = n_cells_x * n_cells_y * n_cells_z; 

	n_faces = (n_cells_x + 1) * n_cells_y *n_cells_z +
		(n_cells_x ) * (n_cells_y + 1) *n_cells_z +
		(n_cells_x ) * n_cells_y *(n_cells_z + 1);

	
	// neighbours is internal faces
	n_neighbours = n_faces - 2 * n_cells_y *n_cells_x
		- 2 * n_cells_x *n_cells_z
		- 2 * n_cells_y *n_cells_z;

	num_bc_cells = n_faces - n_neighbours;

	n_face_x = (n_cells_x - 1) * n_cells_y *n_cells_z;
	n_face_y = (n_cells_x) * (n_cells_y - 1) *n_cells_z;
	n_face_z = (n_cells_x)* n_cells_y *(n_cells_z - 1);

	n_face_points = n_faces * 4;
	owner = new int[n_faces];
	if (owner == NULL) exit(1);
	neighbour = new int[n_faces];
	if (neighbour == NULL) exit(1);


	// get owner neighbour relations and face connectivity

	// start with internal faces
	int face = 0;

	// x-axis faces
	// note n_node = n_cell +1 
	for (int k = 0; k < n_vertices_z -1 ; k++) {
		for (int j = 0; j < n_vertices_y -1;  j++) {
			for (int i = 1; i < n_vertices_x -1 ; i++) {
				
				// first point
				a = k * n_vertices_y* n_vertices_x + j * n_vertices_x  + i;
				b = a + n_vertices_x;
				c = b + n_vertices_y * n_vertices_x;
				d = a + n_vertices_y * n_vertices_x;
				face_connectivity.push_back(std::vector <int>());

				face_connectivity[face].push_back(a);
				face_connectivity[face].push_back(b);
				face_connectivity[face].push_back(c);
				face_connectivity[face].push_back(d);
				owner[face] =  (i-1) + j*n_cells_x + k *n_cells_x *n_cells_y ;
				neighbour[face] = owner[face] +1;

				face++;
			}
		}
	}

	// y- axis faces

	
	for (int k = 0; k < n_vertices_z - 1; k++) {
		for (int j = 1; j < n_vertices_y - 1; j++) {
			for (int i = 0; i < n_vertices_x - 1; i++) {

				// first point
				a = k * n_vertices_y* n_vertices_x + j * n_vertices_x + i;
				b = a + n_vertices_y * n_vertices_x;
				c = b + 1;
				d = a + 1;
				face_connectivity.push_back(std::vector <int>());

				face_connectivity[face].push_back(a);
				face_connectivity[face].push_back(b);
				face_connectivity[face].push_back(c);
				face_connectivity[face].push_back(d);
				owner[face] = (i ) + (j-1) * n_cells_x + k * n_cells_x *n_cells_y;
				neighbour[face] = owner[face] + n_cells_x;

				face++;
			}
		}
	}

	//z-axis faces
	
	for (int k = 1; k < n_vertices_z - 1; k++) {
		for (int j = 0; j < n_vertices_y - 1; j++) {
			for (int i = 0; i < n_vertices_x - 1; i++) {

				// first point
				a = k * n_vertices_y* n_vertices_x + j * n_vertices_x + i;
				b = a + 1;
				c = b + n_vertices_x;
				d = c -1;
				face_connectivity.push_back(std::vector <int>());

				face_connectivity[face].push_back(a);
				face_connectivity[face].push_back(b);
				face_connectivity[face].push_back(c);
				face_connectivity[face].push_back(d);
				owner[face] = (i)+(j ) * n_cells_x + (k-1) * n_cells_x *n_cells_y;
				neighbour[face] = owner[face] + n_cells_x * n_cells_y;

				face++;
			}
		}
	}


	/// boundary faces // left and right

	// x-axis faces
	//left
	for (int k = 0; k < n_vertices_z - 1; k++) {
		for (int j = 0; j < n_vertices_y - 1; j++) {
				i = 0;

				// first point
				a = k * n_vertices_y* n_vertices_x + j * n_vertices_x + i;
				b = a + n_vertices_y * n_vertices_x;
				c = a  + n_vertices_x + n_vertices_y * n_vertices_x;
				d= a + n_vertices_x;

				face_connectivity.push_back(std::vector <int>());

				face_connectivity[face].push_back(a);
				face_connectivity[face].push_back(b);
				face_connectivity[face].push_back(c);
				face_connectivity[face].push_back(d);
				owner[face] = (i ) + j * n_cells_x + k * n_cells_x *n_cells_y;
				neighbour[face] = 0;
				bc_types.push_back("left");
				face++;


			
		}
	}

	for (int k = 0; k < n_vertices_z - 1; k++) {
		for (int j = 0; j < n_vertices_y - 1; j++) {

			i = n_vertices_x -1;
			// first point
			a = k * n_vertices_y* n_vertices_x + j * n_vertices_x + i ;
			b = a + n_vertices_x;
			c = b + n_vertices_y * n_vertices_x;
			d = a + n_vertices_y * n_vertices_x;
			face_connectivity.push_back(std::vector <int>());

			face_connectivity[face].push_back(a);
			face_connectivity[face].push_back(b);
			face_connectivity[face].push_back(c);
			face_connectivity[face].push_back(d);
			owner[face] = (i - 1) + j * n_cells_x + k * n_cells_x *n_cells_y;
			bc_types.push_back("right");
			neighbour[face] = 0;

			face++;


		}
	}

	/// boundary faces // Top and Bottom

	//bottom
	for (int k = 0; k < n_vertices_z - 1; k++) {
	
			for (int i = 0; i < n_vertices_x - 1; i++) {

				j = 0;

				// first point
				a = k * n_vertices_y* n_vertices_x + j * n_vertices_x + i;
				b = a + 1;
				c  = a + n_vertices_y * n_vertices_x + 1; 
				d  = a + n_vertices_y * n_vertices_x;


				face_connectivity.push_back(std::vector <int>());

				face_connectivity[face].push_back(a);
				face_connectivity[face].push_back(b);
				face_connectivity[face].push_back(c);
				face_connectivity[face].push_back(d);
				owner[face] = (i)+(j) * n_cells_x + k * n_cells_x *n_cells_y;
				neighbour[face] = 0;
				bc_types.push_back("bottom");
				face++;
			}
		
	}


	//top
	for (int k = 0; k < n_vertices_z - 1; k++) {
		for (int i = 0; i < n_vertices_x - 1; i++) {

			j = n_vertices_y -1;

			// first point
			a = k * n_vertices_y* n_vertices_x + j * n_vertices_x + i;
			b = a + n_vertices_y * n_vertices_x;
			c = b + 1;
			d = a + 1;
			
			face_connectivity.push_back(std::vector <int>());

			face_connectivity[face].push_back(a);
			face_connectivity[face].push_back(b);
			face_connectivity[face].push_back(c);
			face_connectivity[face].push_back(d);
			owner[face] = (i)+(j - 1) * n_cells_x + k * n_cells_x *n_cells_y;
			neighbour[face] = 0;
			bc_types.push_back("top");
			face++;
		}
	}
	
	///front and back


		//z-axis faces

	
		for (int j = 0; j < n_vertices_y - 1; j++) {
			for (int i = 0; i < n_vertices_x - 1; i++) {

				k = 0;

				// first point
				a = k * n_vertices_y* n_vertices_x + j * n_vertices_x + i;
				b = a + n_vertices_x ; 
				c = b + 1;
				d = a + 1;

			


				
				face_connectivity.push_back(std::vector <int>());

				face_connectivity[face].push_back(a);
				face_connectivity[face].push_back(b);
				face_connectivity[face].push_back(c);
				face_connectivity[face].push_back(d);
				owner[face] = (i)+(j)* n_cells_x + (k ) * n_cells_x *n_cells_y;
				neighbour[face] = 0;
				bc_types.push_back("back");

				face++;
			}
		}

		for (int j = 0; j < n_vertices_y - 1; j++) {
			for (int i = 0; i < n_vertices_x - 1; i++) {

				k = n_vertices_z - 1;

				// first point
				a = k * n_vertices_y* n_vertices_x + j * n_vertices_x + i;
				b = a + 1;
				c = b + n_vertices_x;
				d = c -1;

				


				face_connectivity.push_back(std::vector <int>());

				face_connectivity[face].push_back(a);
				face_connectivity[face].push_back(b);
				face_connectivity[face].push_back(c);
				face_connectivity[face].push_back(d);
				owner[face] = (i)+(j)* n_cells_x + (k-1)* n_cells_x *n_cells_y;
				neighbour[face] = 0;
				bc_types.push_back("front");

				face++;
			}
		}

}



void unstructured_mesh::import_openfoam_mesh(global_variables &globals){
      std::string input_file;


    std::ifstream points;


    std::string line;
    int t,i;

    n_vertices = 0;


    //read in owners file

    input_file = globals.import_file + "owner";

    std::ifstream owner_file;
    owner_file.open(input_file);
    t=0;
    i=0;

    int s1,s2;
    bool start, finish;

    start= false;
    finish = false;

    char comp;

    while (std::getline(owner_file, line))
    {

        std::istringstream iss(line);

        if(t == 12){
            s2 = line.find("nFaces:");
            s1 = line.find("nCells:");

            std::istringstream iss1(line.substr(s1 + 7, s2) );
            iss1 >> n_cells;

            s2 = line.find("nInternalFaces:");
            s1 = line.find("nFaces:");

            std::istringstream iss2(line.substr(s1 + 7, s2) );
            iss2 >> n_faces;

          owner = new int [n_faces];
                if (owner ==NULL) exit (1);
            neighbour = new int [n_faces];
                if (neighbour ==NULL) exit (1);

            std::istringstream iss3(line.substr(s2 + 15, line.length()-1) );
            iss3 >> n_neighbours;



            s2 = line.find("nCells:");
            s1 = line.find("nPoints:");

            std::istringstream iss4(line.substr(s1 + 8, s2) );
            iss4 >> n_vertices;

            x = new float [n_vertices];
                if (x==NULL) exit (1);
            y = new float [n_vertices];
                if (y==NULL) exit (1);
            z = new float [n_vertices];
                if (z==NULL) exit (1);


        }

        if(line.length() >0){
            comp = line.at(0);
        }



        if(comp == ')'){
            start = false;
        }

        if(start && i < n_faces){
            iss >> owner[i] ;
            i++;
        }

        //order of check is important
         if(comp == '('){
            start = true;
        }

        t++;

    }

    owner_file.close();

    input_file = globals.import_file + "points";
    t = 0;
    i =0;
     //read in points file
    points.open(input_file);

    while (std::getline(points, line))
    {

        std::istringstream iss(line);

        if(line.length() >0){
            comp = line.at(0);
        }

        if(comp == ')'){
            start = false;
        }

        if(start){
            iss.ignore(1);
            iss >> x[i] >> y[i] >> z[i] ;
            openfoam_map.push_back(i);
            i++;
        }

        //order of check is important
         if(comp == '('){
            start = true;
        }

        t++;
        // process pair (a,b)
    }

    points.close();



        //read in neighbours file

    input_file = globals.import_file + "neighbour";

    std::ifstream neighbour_file;
    neighbour_file.open(input_file);
    t=0;
    i=0;


    while (std::getline(neighbour_file, line))
    {

        std::istringstream iss(line);

         if(line.length() >0){
            comp = line.at(0);
        }


        if(comp == ')'){
            start = false;
        }

        if(start){
            iss >> neighbour[i] ;
            i++;
        }

        if(comp == '('){
            start = true;
        }
        t++;

    }

    neighbour_file.close();

    for(int i = n_neighbours; i< n_faces;i++){
        neighbour[i] =0;

    }

         //read in neighbours file

    input_file = globals.import_file + "faces";

    std::ifstream faces_file;
    faces_file.open(input_file);
    t=0;
    i=0;

    n_face_points =0;
    int a, b,c,d,e,f;
    while (std::getline(faces_file, line))
    {

        std::istringstream iss(line);


        if(t == 17){
            iss >> n_faces;
         }


      if(line.length() >0){
            comp = line.at(0);
        }


        if(comp == ')'){
            start = false;
        }

        if(start){

            iss >> a;

            n_face_points = n_face_points + a;
            iss.ignore(1);

            face_connectivity.push_back(std::vector <int> ());
			if (a == 3) {
				iss >> a >> b >> c ;
				face_connectivity[i].push_back(a);
				face_connectivity[i].push_back(b);
				face_connectivity[i].push_back(c);
			}else if (a == 4) {
                 iss >> a >> b >> c >> d ;
                 face_connectivity[i].push_back(a);
                face_connectivity[i].push_back(b);
                face_connectivity[i].push_back(c);
                face_connectivity[i].push_back(d);

            }else if (a == 5){
                 iss >> a >> b >> c >> d >> e;
                face_connectivity[i].push_back(a);
                face_connectivity[i].push_back(b);
                face_connectivity[i].push_back(c);
                face_connectivity[i].push_back(d);
                face_connectivity[i].push_back(e);
            }else if (a == 6){
                iss >> a >> b >> c >> d >> e >>f;
                  face_connectivity[i].push_back(a);
                face_connectivity[i].push_back(b);
                face_connectivity[i].push_back(c);
                face_connectivity[i].push_back(d);
                face_connectivity[i].push_back(e);
                face_connectivity[i].push_back(f);
            }
            i++;
        }

         if(comp == '('){
            start = true;
        }


        t++;



        t++;

    }

    faces_file.close();


    num_bc_cells = n_faces - n_neighbours;


}


void unstructured_mesh::output_mesh_to_text(global_variables &globals){

    std::ofstream mesh_output ;
    std::string output_dir;

    output_dir = globals.output_file +"/mesh.txt";


        // error_output.open("/home/brendan/Dropbox/PhD/Test Cases/Couette Flow/error.txt", ios::out);
    mesh_output.open(output_dir.c_str(), ios::out);
    for (int mesh_t= 0; mesh_t < total_cells; mesh_t++){
         mesh_output << mesh_t << "," << centroid_x[mesh_t] << "," <<
                centroid_y[mesh_t] << "," << centroid_z[mesh_t] << "," << cell_volume[mesh_t]

                 << endl;


        }

     for (int mesh_t= 0; mesh_t < n_faces; mesh_t++){
         mesh_output << mesh_t << "," << face_area[mesh_t] <<
                "," << face_x[mesh_t] << "," << face_y[mesh_t] << "," << face_z[mesh_t] <<
                 "," << face_i[mesh_t] << "," << face_j[mesh_t] << "," << face_k[mesh_t] << ","  << delta_t_face[mesh_t] << ","
                 << endl;


        }


    mesh_output.close();
}

void unstructured_mesh::calc_centre_node(){

    double c_x = X/2;
    double c_y = Y/2;
    double c_z = Z/2;


    double min_dist,dist;

    min_dist = X*Y*Z;

    vector_var domain_centre, cell;
    domain_centre.set_equal(c_x,c_y,c_z);
    centre_node = 0;

    for(int i = 0; i< n_cells;i++){

        cell.set_equal(centroid_x[i],centroid_y[i],centroid_z[i]);

        cell.subtract(domain_centre);

        dist = cell.Magnitude();

        if (dist < min_dist){
            min_dist = dist;
            centre_node = i;
        }




    }

    return;
}




void unstructured_mesh::tet_volume (int q, int r, int s, std::vector<vector_var> &nodes, int i){

    double vol;

    //http://paulbourke.net/geometry/polygonmesh/

    vol =0.0;

    vector_var a,b,c;

    // copy vctor_vars from vector
    //negatives are to allow for zero indexing
    a = nodes[q];
    b = nodes[r];
    c = nodes[s];

    b.subtract(a);
    c.subtract(a);
    //cross product
    b.cross_product(c);
    //dot productto find volume
    vol = a.Dot_Product(b)/6;
    // add contribution to cell volume
    cell_volume[owner[i]] = cell_volume[owner[i]] + vol;

    if(i < n_neighbours ){
        cell_volume[neighbour[i]] = cell_volume[neighbour[i]] - vol;
    }


}

//REDUNDANT FUNCTION

//void unstructured_mesh::calc_cell_volume(std::vector<vector_var> &nodes, int i){
//
//
//
//    cell_volume[i] = 0;
//    //split into 6 pyramids
//
//    //node notations, 3 square pyramids split into 2
//
//    //notation is CGNS HEXA_8
//
//    // 71256
//    // 7 125
//    //tet_volume( 1,2,5,7,nodes,i);
//    // 7 526
//    //tet_volume( 6,2,5,7,nodes,i);
//
//    //78541
//    // 7 854
//    tet_volume( 8,5,4,7,nodes,i);
//    // 7541
//    tet_volume( 5,4,1,7,nodes,i);
//
//    //71234
//    // 7 123
//    tet_volume( 1,2,3,7,nodes,i);
//    // 7 134
//    tet_volume( 1,3,4,7,nodes,i);
//
//}

void unstructured_mesh::remove_face_from_flux_calcs(int i, int n){

        if (w_node[n] ==i){
                calculated_faces[n].erase(
                    std::remove(calculated_faces[n].begin(), calculated_faces[n].end(), 0)
                    , calculated_faces[n].end());
             }else if (s_node[n] ==i){
                calculated_faces[n].erase(
                    std::remove(calculated_faces[n].begin(), calculated_faces[n].end(), 1)
                    , calculated_faces[n].end());

            }else if (b_node[n] ==i){
                calculated_faces[n].erase(
                    std::remove(calculated_faces[n].begin(), calculated_faces[n].end(), 5)
                    , calculated_faces[n].end());

            }else if (n_node[n] ==i){
                calculated_faces[n].erase(
                    std::remove(calculated_faces[n].begin(), calculated_faces[n].end(), 3)
                    , calculated_faces[n].end());

            }else if (e_node[n] ==i){
                calculated_faces[n].erase(
                    std::remove(calculated_faces[n].begin(), calculated_faces[n].end(), 2)
                    , calculated_faces[n].end());
            }else if (f_node[n] ==i){
                calculated_faces[n].erase(
                    std::remove(calculated_faces[n].begin(), calculated_faces[n].end(), 4)
                    , calculated_faces[n].end());

            }



}


void unstructured_mesh::calculate_face_flux_calcs(){


    // initialise vector of vectors
    // interior cells should orinally have 6 calcs and calcs in all 6 directions
    for( int i = 0; i < n_cells; i++){

        calcs_per_cell[i] = 6;
        calculated_faces.push_back(std::vector <int> ());
        calculated_faces[i].push_back(0);//west
        calculated_faces[i].push_back(1); // south
        calculated_faces[i].push_back(2); // east
        calculated_faces[i].push_back(3);// north
        calculated_faces[i].push_back(4); // front
        calculated_faces[i].push_back(5); // back

    }


    //remove duplicate flux calcs

    int n,e,f;

    for( int i = 0; i < n_cells; i++){

            n = n_node[i];
            if( n >-1){
                // reduce calcs by 1 in neighbouring cell
                calcs_per_cell[n] = calcs_per_cell[n] -1;
                remove_face_from_flux_calcs(i,n);
            }

            e = e_node[i];
            if( e >-1){
                // reduce calcs by 1 in neighbouring cell
                calcs_per_cell[e] = calcs_per_cell[e] -1;
                remove_face_from_flux_calcs(i,e);
            }

            f = f_node[i];
            if( f >-1){
                // reduce calcs by 1 in neighbouring cell
                calcs_per_cell[f] = calcs_per_cell[f] -1;
                remove_face_from_flux_calcs(i,f);
            }

    }


}



void unstructured_mesh::cell_tet_area(std::vector<vector_var> &nodes, int i, int q, int r, int s,
        double &face_area, double &total_area,vector_var &centroid,vector_var &face_centroid,
        vector_var &face_normal){

        vector_var a,b,c,d,normal;
        double area;
        //http://paulbourke.net/geometry/polygonmesh/
        // note that area calcs here are missing 0.5 factor as they cancel in overall calc

        a = nodes[q];
        b = nodes[r];
        c = nodes[s];

        //centroid of tet
        d.average(a,b,c);
        b.subtract(a);
        c.subtract(a);
        b.cross_product(c);
        area =  0.5* b.Magnitude();
        face_area = face_area + area;
        total_area = total_area + area;

        d.factor(1/area);
        centroid.add(d);
        face_centroid.add(d);

        face_normal.add(b) ;




        //facenormal check for inward/outward


}


void unstructured_mesh::calc_face_centroid(std::vector<vector_var> &nodes, int i){

    vector_var a,b,c,d, centroid,face_centroid;
    double total_area;

}

void unstructured_mesh::calc_centroids_volume_area_normal(std::vector<vector_var> &nodes, int i){

    vector_var centroid,face_centroid, face_normal;
    vector_var zero(0.0,0.0,0.0);

    double total_area;

   n_area[i] =0.0;
   e_area[i] =0.0;
   w_area[i] =0.0;
   s_area[i] =0.0;
   f_area[i] =0.0;
   b_area[i] =0.0;

    total_area = 0.0;
    //http://paulbourke.net/geometry/polygonmesh/

    /// note order of arguments of tets gives correct outward normal
   ///north area
    //3478
    //split into 2 pyramids
    face_centroid.set_equal(zero);
    face_normal.set_equal(zero);
    centroid.set_equal(zero);

    cell_tet_area(nodes,i, 7,4,3,n_area[i],total_area,centroid,face_centroid,face_normal);
    cell_tet_area(nodes,i, 7,4,8,n_area[i],total_area,centroid,face_centroid,face_normal);
    north_x[i] = face_centroid.x/n_area[i];
    north_y[i] = face_centroid.y/n_area[i];
    north_z[i] = face_centroid.z/n_area[i];

    n_i[i] = face_normal.x;
    n_j[i] =face_normal.y;
    n_k[i] = face_normal.z;

///south area
    //1256
    //split into 2 pyramids
    face_centroid.set_equal(zero);
     cell_tet_area(nodes,i, 1,6,2,s_area[i],total_area,centroid,face_centroid,face_normal);
    cell_tet_area(nodes,i, 1,6,5,s_area[i],total_area,centroid,face_centroid,face_normal);
    south_x[i] = face_centroid.x/s_area[i];
    south_y[i] = face_centroid.y/s_area[i];
    south_z[i] = face_centroid.z/s_area[i];
    s_i[i] = face_normal.x;
    s_j[i] =face_normal.y;
    s_k[i] = face_normal.z;

///east area
    //2376
    //split into 2 pyramids
    face_centroid.set_equal(zero);
    cell_tet_area(nodes,i, 6,7,2,e_area[i],total_area,centroid,face_centroid,face_normal);
    cell_tet_area(nodes,i, 3,7,2,e_area[i],total_area,centroid,face_centroid,face_normal);
    east_x[i] = face_centroid.x/e_area[i];
    east_y[i] = face_centroid.y/e_area[i];
    east_z[i] = face_centroid.z/e_area[i];
    e_i[i] = face_normal.x;
    e_j[i] =face_normal.y;
    e_k[i] = face_normal.z;

///west area
    //1584
    //split into 2 pyramids
    face_centroid.set_equal(zero);
    cell_tet_area(nodes,i, 4,5,8,w_area[i],total_area,centroid,face_centroid,face_normal);
    cell_tet_area(nodes,i, 1,5,4,w_area[i],total_area,centroid,face_centroid,face_normal);
    west_x[i] = face_centroid.x/w_area[i];
    west_y[i] = face_centroid.y/w_area[i];
    west_z[i] = face_centroid.z/w_area[i];
    w_i[i] = face_normal.x;
    w_j[i] =face_normal.y;
    w_k[i] = face_normal.z;

    ///front area
    //5678
    //split into 2 pyramids
    face_centroid.set_equal(zero);
    cell_tet_area(nodes,i, 6,5,7,f_area[i],total_area,centroid,face_centroid,face_normal);
    cell_tet_area(nodes,i, 8,5,7,f_area[i],total_area,centroid,face_centroid,face_normal);
    front_x[i] = face_centroid.x/f_area[i];
    front_y[i] = face_centroid.y/f_area[i];
    front_z[i] = face_centroid.z/f_area[i];
    f_i[i] = face_normal.x;
    f_j[i] =face_normal.y;
    f_k[i] = face_normal.z;

///back area
    //1234
    //split into 2 pyramids
    face_centroid.set_equal(zero);
    cell_tet_area(nodes,i, 2,3,1,b_area[i],total_area,centroid,face_centroid,face_normal);
    cell_tet_area(nodes,i, 4,3,1,b_area[i],total_area,centroid,face_centroid,face_normal);
    back_x[i] = face_centroid.x/b_area[i];
    back_y[i] = face_centroid.y/b_area[i];
    back_z[i] = face_centroid.z/b_area[i];
    b_i[i] = face_normal.x;
    b_j[i] =face_normal.y;
    b_k[i] = face_normal.z;


    centroid.factor(total_area);

    centroid_x[i] = centroid.x;
    centroid_y[i] = centroid.y;
    centroid_z[i] = centroid.z;

}

void unstructured_mesh::calc_internal_sphere_radius(){

    int cell_a, cell_b;
    for(int i = 0; i < n_faces; i++){

        // get distance from face_centroid to cell_centroid

        cell_a = owner[i];
        cell_b = neighbour[i];


        vector_var face(face_x[i],face_y[i],face_z[i]);

        vector_var centroid(centroid_x[cell_a],centroid_y[cell_a],centroid_z[cell_a]);
        centroid.subtract(face);
        cell_cfl_r[cell_a] = min(cell_cfl_r[cell_a] ,centroid.Magnitude());

          centroid.set_equal(centroid_x[cell_b],centroid_y[cell_b],centroid_z[cell_b]);
        centroid.subtract(face);
        cell_cfl_r[cell_b] = min(cell_cfl_r[cell_b] ,centroid.Magnitude());

    }


}

void unstructured_mesh::calc_face_delta_t(std::vector<vector_var> &nodes, int i,global_variables &globals){

    //algorithm is to find min of delta _x, y _z
    // the find skewness
    // delta t = 0.5 * min_delta * cos(skewness)

    double angle, min_angle, max_angle, skewness,min_delta;

    //angles in radians
    min_angle = 2;
    max_angle = 0;



     double max_delta_x,max_delta_y,max_delta_z;


    //min dx

    vector<double> nodes_x ;

    for ( int j =0; j < nodes.size(); j++){
        nodes_x.push_back(nodes[j].x);
    }
    double dx1,dx2;

    dx1 = fabs(nodes_x[0] - nodes_x[2]);
    dx2 = fabs(nodes_x[0] - nodes_x[2]);
    max_delta_x =  max(dx1,dx2);
    //find min-non zero value; greater than 0.1 is limit

    vector<double> nodes_y;
    for ( int j =0; j < nodes.size(); j++){
        nodes_y.push_back(nodes[j].y);
    }

     double dy1,dy2;

    dy1 = fabs(nodes_y[0] - nodes_y[2]);
    dy2 = fabs(nodes_y[0] - nodes_y[2]);
    max_delta_y =  max(dy1,dy2);



    vector<double> nodes_z;
    for ( int j =0; j < nodes.size(); j++){
            nodes_z.push_back(nodes[j].z);
        }

    double dz1,dz2;

    dz1 = fabs(nodes_z[0] - nodes_z[2]);
    dz2 = fabs(nodes_z[0] - nodes_z[2]);
    max_delta_z =  max(dz1,dz2);

    cell_cfl_x[owner[i]] = max(cell_cfl_x[owner[i]], max_delta_x);
    cell_cfl_x[neighbour[i]] = max(cell_cfl_x[owner[i]], max_delta_x);

    cell_cfl_y[owner[i]] = max(cell_cfl_y[owner[i]], max_delta_y);
    cell_cfl_y[neighbour[i]] = max(cell_cfl_y[owner[i]], max_delta_y);

    cell_cfl_z[owner[i]] = max(cell_cfl_z[owner[i]], max_delta_z);
    cell_cfl_z[neighbour[i]] = max(cell_cfl_z[owner[i]], max_delta_z);


	double reference_dt, face_dt;

	// need to get accurate dt
	// check for distances to centroid of adjoining volumes
	// then check for square root of area
	// then stability requirements

	/*vector_var face(face_x[i], face_y[i], face_z[i]);

	vector_var centroid(centroid_x[owner[i]], centroid_y[owner[i]], centroid_z[owner[i]]);
	centroid.subtract(face);
	face_dt = centroid.Magnitude();


	centroid.set_equal(centroid_x[neighbour[i]], centroid_y[neighbour[i]], centroid_z[neighbour[i]]);
	centroid.subtract(face);

	face_dt = min(face_dt, centroid.Magnitude());

	face_dt = min(face_dt, 0.5 * sqrt(face_area[i]));*/

	face_dt = 0.5* sqrt(face_area[i]);

	// minimum tau for stability is 0.55 (conservative) or user input
	 reference_dt = globals.visc * 3 / globals.pre_conditioned_gamma / globals.fneq_min;
	 
	 //for inviscid flow, relation doesn't hold
	 if (globals.testcase != 8) {
		 delta_t_face[i] = min(face_dt, reference_dt);
	 }
    i =i +1;

}



void unstructured_mesh::calc_delta_t(std::vector<vector_var> &nodes, int i,global_variables &globals){

    //algorithm is to find min of delta _x, y _z
    // the find skewness
    // delta t = 0.5 * min_delta * cos(skewness)
    vector_var n (n_i[i], n_j[i], n_k[i]);
    vector_var s (s_i[i], s_j[i], s_k[i]);
    vector_var e (e_i[i], e_j[i], e_k[i]);
    vector_var w (w_i[i], w_j[i], w_k[i]);
    vector_var f (f_i[i], f_j[i], f_k[i]);
    vector_var b (b_i[i], b_j[i], b_k[i]);

    double angle, min_angle, max_angle, skewness;

    //angles in radians
    min_angle = 2;
    max_angle = 0;

    // 6 planes in hex, lead to 12 dihedral angles

    //ne
    angle = acos(n.Dot_Product(e)/n.Magnitude()/e.Magnitude());
    min_angle = min(angle,min_angle);
    max_angle = max(angle, max_angle);
       //nw
    angle = acos(n.Dot_Product(w)/n.Magnitude()/w.Magnitude());
    min_angle = min(angle,min_angle);
    max_angle = max(angle, max_angle);
    //nf
    angle = acos(n.Dot_Product(f)/n.Magnitude()/f.Magnitude());
    min_angle = min(angle,min_angle);
    max_angle = max(angle, max_angle);
    //nb
    angle = acos(n.Dot_Product(b)/n.Magnitude()/b.Magnitude());
    min_angle = min(angle,min_angle);
    max_angle = max(angle, max_angle);

    //se
     angle = acos(s.Dot_Product(e) /s.Magnitude()/e.Magnitude());
    min_angle = min(angle,min_angle);
    max_angle = max(angle, max_angle);

    //sw

     angle = acos(s.Dot_Product(w)/s.Magnitude()/w.Magnitude());
    min_angle = min(angle,min_angle);
    max_angle = max(angle, max_angle);
    //sf

     angle = acos(s.Dot_Product(f) /s.Magnitude()/f.Magnitude());
    min_angle = min(angle,min_angle);
    max_angle = max(angle, max_angle);
    //sb

     angle = acos(s.Dot_Product(b) /s.Magnitude()/b.Magnitude());
    min_angle = min(angle,min_angle);
    max_angle = max(angle, max_angle);
    //fe

     angle = acos(f.Dot_Product(e) /f.Magnitude()/e.Magnitude());
    min_angle = min(angle,min_angle);
    max_angle = max(angle, max_angle);
    //fw

     angle = acos(f.Dot_Product(w) /f.Magnitude()/w.Magnitude());
    min_angle = min(angle,min_angle);
    max_angle = max(angle, max_angle);
    //be

     angle = acos(b.Dot_Product(e) /b.Magnitude()/e.Magnitude());
    min_angle = min(angle,min_angle);
    max_angle = max(angle, max_angle);
    //bw

     angle = acos(b.Dot_Product(w) /b.Magnitude()/w.Magnitude());
    min_angle = min(angle,min_angle);
    max_angle = max(angle, max_angle);

    //90 degrees in ideal dihedral angle for hexahedron = pi/2
    skewness = max ( ( max_angle - globals.PI/2)/globals.PI*2 , (1 - min_angle/globals.PI*2) );

    double min_delta;
    min_delta = 0.0;

    //min dx

    vector<double> nodes_x ;

    for ( int j =0; j < nodes.size(); j++){
        nodes_x.push_back(nodes[j].x);
    }

    sort(nodes_x.begin(), nodes_x.end());
    adjacent_difference(begin(nodes_x), end(nodes_x), begin( nodes_x));
    sort(nodes_x.begin(), nodes_x.end());
    auto lower = lower_bound(nodes_x.begin(),nodes_x.end(), 0.01);
    //find min-non zero value; greater than 0.1 is limit
    min_delta =  nodes_x[distance(nodes_x.begin(),lower)];


    vector<double> nodes_y;
    for ( int j =0; j < nodes.size(); j++){
        nodes_y.push_back(nodes[j].y);
    }

    sort(nodes_y.begin(), nodes_y.end());
    adjacent_difference(begin(nodes_y), end(nodes_y), begin(nodes_y));
    sort(nodes_y.begin(), nodes_y.end());

 lower = lower_bound(nodes_y.begin(),nodes_y.end(), 0.01);
    //find min-non zero value; greater than 0.1 is limit
    min_delta = min(min_delta, nodes_y[distance(nodes_y.begin(),lower)] );



    vector<double> nodes_z;
    for ( int j =0; j < nodes.size(); j++){
            nodes_z.push_back(nodes[j].z);
        }

    sort(nodes_z.begin(), nodes_z.end());
    adjacent_difference(begin(nodes_z), end(nodes_z), begin(nodes_z));
    sort(nodes_z.begin(), nodes_z.end());

    lower = lower_bound(nodes_z.begin(),nodes_z.end(), 0.01);
    //find min-non zero value; greater than 0.1 is limit
    min_delta = min(min_delta,nodes_z[distance(nodes_z.begin(),lower)]);

    min_delta = min( e_area[i], w_area[i]);
    min_delta = min( min_delta, n_area[i]);
    min_delta = min( min_delta, s_area[i]);
    min_delta = min( min_delta, f_area[i]);
    min_delta = min( min_delta, b_area[i]);

    delta_t[i] = 0.5 * sqrt(min_delta) * cos(skewness);
    i =i +1;

}




void unstructured_mesh::generate_internal_mesh(global_variables &globals){


    ///get volume contribution of each tet
    std::vector <vector_var> nodes;
    nodes.push_back(vector_var(0,0,0));
    int t,nb;
    double total_area,vol;
    vector_var face_centroid, face_normal,centroid;

    for(int i = 0; i< total_cells; i++){

        gradient_cells.push_back(std::vector <int> ());
        gradient_faces.push_back(std::vector <int> ());
    }

     for(int i =0 ; i< n_faces ; i ++){
        nodes.clear();
        for(int j = 0; j<face_connectivity[i].size(); j ++){
            t =face_connectivity[i][j];
            nodes.push_back(vector_var(x[t],y[t],z[t]));
            //nodes.push_back(vector_var(round(x[t]),round(y[t]),round(z[t])));
        }

        for( int t = 0; t< face_connectivity[i].size() -2; t++){
         tet_volume(0,(t+1),(t+2),nodes,i);
        }

        total_area = 0.0;
        centroid.set_equal(0.0,0.0,0.0);
        face_centroid.set_equal(0.0,0.0,0.0);
        face_normal.set_equal(0.0,0.0,0.0);

        for( int t = 0; t< face_connectivity[i].size() -2; t++){
            cell_tet_area(nodes,i, 0,(t+1),(t+2),face_area[i],total_area,centroid,face_centroid,face_normal);
        }

        total_cell_area[owner[i]] = total_cell_area[owner[i]]  + face_area[i];

        centroid_x[owner[i]] = centroid_x[owner[i]] + face_centroid.x;
        centroid_y[owner[i]] = centroid_y[owner[i]] + face_centroid.y;
        centroid_z[owner[i]] = centroid_z[owner[i]] + face_centroid.z;

        face_x[i] = face_centroid.x/face_area[i];
        face_y[i] = face_centroid.y/face_area[i];
        face_z[i] = face_centroid.z/face_area[i];

        face_i[i] = face_normal.x/face_normal.Magnitude();
        face_j[i] = face_normal.y/face_normal.Magnitude();
        face_k[i] = face_normal.z/face_normal.Magnitude();

        calc_face_delta_t(nodes,i,globals);

        if( i > (n_neighbours -1)) {

            nb = i - n_neighbours +n_cells ;
            gradient_cells[owner[i]].push_back(nb);
            gradient_faces[owner[i]].push_back(i);

        }else{

            total_cell_area[neighbour[i]] = total_cell_area[neighbour[i]]  + face_area[i];

            centroid_x[neighbour[i]] = centroid_x[neighbour[i]] + face_centroid.x;
            centroid_y[neighbour[i]] = centroid_y[neighbour[i]] + face_centroid.y;
            centroid_z[neighbour[i]] = centroid_z[neighbour[i]] + face_centroid.z;

            gradient_cells[owner[i]].push_back(neighbour[i]);
            gradient_cells[neighbour[i]].push_back(owner[i]);

            gradient_faces[owner[i]].push_back(i);
            gradient_faces[neighbour[i]].push_back(i);

        }


    }

    for(int i=0; i< n_cells; i++){

        centroid_x[i] = centroid_x[i]  / total_cell_area[i];
        centroid_y[i] = centroid_y[i]  / total_cell_area[i];
        centroid_z[i] = centroid_z[i]  / total_cell_area[i];

    }

    //get centroids of ghost cells
    for(int i = n_neighbours; i < n_faces; i++){
        populate_ghost_cell_from_foam_face(i);

    }


    calc_internal_sphere_radius();

    calc_centre_node();



    return;
}
void unstructured_mesh::destruct_mesh_variables(global_variables &globals){

       if(globals.import_format == "openfoam"){

        }else{
//            delete[] ielem[0];
//            delete[] ielem;
//            delete[] ielem_mix[0];
//            delete[] ielem_mix;
//
//            delete[] i_elem_bc[0];
//            delete[] i_elem_bc;
//            delete[] iparentdata[0];
//            delete[] iparentdata;

        }

	   std::cout << "Destructing Intermediate Variables" << endl;

	   delete[](x);
	   x = NULL;
	   delete[](y);
	   y = NULL;
	   delete[](z);
	   z = NULL;

	   delete[](front_x);
	   front_x = NULL;
	   delete[](front_y);
	   front_y = NULL;
	   delete[](front_z);
	   front_z = NULL;
	   delete[](back_x);
	   back_x = NULL;
	   delete[](back_y);
	   back_y = NULL;
	   delete[](back_z);
	   back_z = NULL;


	   delete[](f_area);
	   f_area = NULL;
	   delete[](f_i);
	   f_i = NULL;
	   delete[](f_j);
	   f_j = NULL;
	   delete[](f_k);
	   f_k = NULL;
	   delete[](f_node);
	   f_node = NULL;

	   delete[](b_area);
	   b_area = NULL;
	   delete[](b_i);
	   b_i = NULL;
	   delete[](b_j);
	   b_j = NULL;
	   delete[](b_k);
	   b_k = NULL;
	   delete[](b_node);
	   b_node = NULL;

	   delete[](cell_cfl_x);
	   cell_cfl_x = NULL;

	   delete[](cell_cfl_y);
	   cell_cfl_y = NULL;

	   delete[](cell_cfl_z);
	   cell_cfl_z = NULL;

	   delete[](cell_cfl_r);
	   cell_cfl_r = NULL;

	   //dtor
	 
	   delete[](north_x);
	   north_x = NULL;
	   delete[](south_x);
	   south_x = NULL;
	   delete[](east_x);
	   east_x = NULL;
	   delete[](west_x);
	   west_x = NULL;


	   
	   delete[](north_y);
	   north_y = NULL;
	   delete[](south_y);
	   south_y = NULL;
	   delete[](east_y);
	   east_y = NULL;
	   delete[](west_y);
	   west_y = NULL;

	  
	   delete[](north_z);
	   north_z = NULL;
	   delete[](south_z);
	   south_z = NULL;
	   delete[](east_z);
	   east_z = NULL;
	   delete[](west_z);
	   west_z = NULL;

	  
	  
	   delete[](delta_t_n);
	   delta_t_n = NULL;
	   delete[](delta_t_e);
	   delta_t_e = NULL;

	   delete[](n_area);
	   n_area = NULL;
	   delete[](n_i);
	   n_i = NULL;
	   delete[](n_j);
	   n_j = NULL;
	   delete[](n_k);
	   n_k = NULL;
	   delete[](n_node);
	   n_node = NULL;

	   delete[](e_area);
	   e_area = NULL;
	   delete[](e_i);
	   e_i = NULL;
	   delete[](e_j);
	   e_j = NULL;
	   delete[](e_k);
	   e_k = NULL;
	   delete[](e_node);
	   e_node = NULL;

	   delete[](w_area);
	   w_area = NULL;
	   delete[](w_i);
	   w_i = NULL;
	   delete[](w_j);
	   w_j = NULL;
	   delete[](w_k);
	   w_k = NULL;
	   delete[](w_node);
	   w_node = NULL;

	   delete[](s_area);
	   s_area = NULL;
	   delete[](s_i);
	   s_i = NULL;
	   delete[](s_j);
	   s_j = NULL;
	   delete[](s_k);
	   s_k = NULL;
	   delete[](s_node);
	   s_node = NULL;


	      calculated_faces.clear();
	   calculated_faces.shrink_to_fit();

	   face_connectivity.clear();
	   face_connectivity.shrink_to_fit();

	   volume_connectivity.clear();
	   volume_connectivity.shrink_to_fit();

	   ghost_faces.clear();
	   ghost_faces.shrink_to_fit();

	   face_labels.clear();
	   face_labels.shrink_to_fit();

	   std::cout << "Finished Destructing Variables" << endl;


}


void unstructured_mesh::calc_unstructured_lattice_size(){


    for (int i = 0; i< n_cells; i++){
            //get the north and east flux streaming time step
            if (n_node[i] > -1 ){
                delta_t_n[i] = min(delta_t[i], delta_t[n_node[i]]);
            }
            if (e_node[i] > -1 ){
                delta_t_e[i] = min(delta_t[i], delta_t[e_node[i]]);
            }
    }

}

void unstructured_mesh::get_ghost_face_hash(global_variables &globals){

int neighbour;

std::vector <int> temp;
long long int hash_int;
int t,of_face;
int factor;

factor = ceil(log10(n_cells));

for (int i = 0; i < num_bc_cells; i++){

    temp.clear();


    if(globals.import_format == "openfoam"){
        of_face = i + n_neighbours;
        neighbour = owner[of_face];

        temp.push_back(face_connectivity[of_face][0]);
        temp.push_back(face_connectivity[of_face][1]);
        temp.push_back(face_connectivity[of_face][2]);
        temp.push_back(face_connectivity[of_face][3]);


    }else{
        // need neighbouring cells
//        neighbour = iparentdata[0][i];
//
//        temp.push_back(i_elem_bc[0][i*5 +1]);
//        temp.push_back(i_elem_bc[0][i*5 +2]);
//        temp.push_back(i_elem_bc[0][i*5 +3]);
//        temp.push_back(i_elem_bc[0][i*5 +4]);


    }

    std::sort( temp.begin(), temp.end());
    hash_int =0;
    t = 0;

    for(std::vector<int>::iterator it = temp.begin(); it != temp.end(); ++it) {

       hash_int = hash_int + *it * pow(10,factor * t);
       t++;
    }

    ghost_faces.push_back(hash_int);
    // then identify the ghost cell parameters needed

}

}

void unstructured_mesh::populate_ghost_cells(global_variables &globals){


    get_ghost_face_hash(globals);
    std::vector<long long int>::iterator it;
    int p;
    //find position in face
     for(int i =0; i < num_bc_cells; i++){
            it = std::find(face_labels.begin(),face_labels.end(),ghost_faces[i]);
            if( it!= face_labels.end()){
                p = std::distance( face_labels.begin(), it );
                populate_ghost_face_neighbour(i,p);
            }



     }
    // get diatance

}

void unstructured_mesh::populate_ghost_cell_from_foam_face(int face){

    int ghost_element,element2,k,j,nb;

    nb = n_cells + face- n_neighbours;

    vector_var face_centroid ;

    face_centroid.x =face_x[face];
    face_centroid.y =face_y[face];
    face_centroid.z = face_z[face];

    ///populate

    j = owner[face];

    vector_var neighbour_centroid(centroid_x[j], centroid_y[j], centroid_z[j]);

    //get relative vector
    face_centroid.subtract(neighbour_centroid);
    //multiply by 2
    face_centroid.factor(0.5);

    /*centroid_x[nb] = neighbour_centroid.x + face_centroid.x * fabs(face_i[face]) ;
    centroid_y[nb] = neighbour_centroid.y + face_centroid.y * fabs(face_j[face]);
    centroid_z[nb] = neighbour_centroid.z + face_centroid.z * fabs(face_k[face]);*/

	centroid_x[nb] = neighbour_centroid.x + face_centroid.x ;
	centroid_y[nb] = neighbour_centroid.y + face_centroid.y ;
	centroid_z[nb] = neighbour_centroid.z + face_centroid.z ;

    //add bc neighbour
    neighbour[face]  = nb;
}


void unstructured_mesh::populate_ghost_face_neighbour(int ghost, int index2){

    int ghost_element,element2,k;

    ghost_element = ghost + n_cells ;
    element2 = floor(index2/6);

    int face_i[6] = {1,0,1,0,0,0};
    int face_j[6] = {0,1,0,1,0,0};
    int face_k[6] = {0,0,0,0,1,1};

    k = index2%6;
    vector_var face_centroid ;
    switch(k) {
            case 0: // west
                w_node[element2] = ghost_element;

                face_centroid.x = west_x[element2];
                 face_centroid.y = west_y[element2];
                 face_centroid.z = west_z[element2];

                 break;
            case 1: // south
                s_node[element2] = ghost_element;

                face_centroid.x = south_x[element2];
                face_centroid.y = south_y[element2];
                face_centroid.z = south_z[element2];
                 break;
            case 2: // east
                e_node[element2] = ghost_element;
                face_centroid.x = east_x[element2];
                face_centroid.y = east_y[element2];
                face_centroid.z = east_z[element2];
                 break;
            case 3: // north
                n_node[element2] = ghost_element;
                face_centroid.x = north_x[element2];
                face_centroid.y = north_y[element2];
                face_centroid.z = north_z[element2];
                 break;
            case 4: // front
                f_node[element2] = ghost_element;
                face_centroid.x =front_x[element2];
                face_centroid.y = front_y[element2];
                face_centroid.z = front_z[element2];
                 break;
            case 5: // back
                b_node[element2] = ghost_element;
                face_centroid.x =back_x[element2];
                face_centroid.y =back_y[element2];
                face_centroid.z = back_z[element2];
                 break;
        }

    ///populate

    bc_face_x[ghost] = face_centroid.x;
    bc_face_y[ghost] = face_centroid.y;
    bc_face_z[ghost] = face_centroid.z;

    bc_face_direction[ghost] = k;

    vector_var neighbour_centroid(centroid_x[element2], centroid_y[element2], centroid_z[element2]);

    //get relative vector
    face_centroid.subtract(neighbour_centroid);
    //multiply by 2
    face_centroid.factor(0.5);

    centroid_x[ghost_element] = neighbour_centroid.x +face_centroid.x * face_i[k] ;
    centroid_y[ghost_element] = neighbour_centroid.y + face_centroid.y * face_j[k];
    centroid_z[ghost_element] = neighbour_centroid.z + face_centroid.z * face_k[k];

    delta_t[ghost_element] = delta_t[element2];

    //add bc neighbour
    bc_neighbour[ghost]  = element2;
}


void unstructured_mesh::face_neighbours_hash(int a, int b, int c, int d, int i, global_variables &globals){

std::vector <int> temp;
int ijk_to_cgns[8] = {0,4, 6,2,1,5,7,3};

if(globals.import_format == "openfoam"){

    temp.push_back(volume_connectivity[i][ijk_to_cgns[a-1]]);
    temp.push_back(volume_connectivity[i][ijk_to_cgns[b-1]]);
    temp.push_back(volume_connectivity[i][ijk_to_cgns[c-1]]);
    temp.push_back(volume_connectivity[i][ijk_to_cgns[d-1]]);

}else{
    //get 4 vertices into vector
//    temp.push_back(ielem_mix[0][i*9 + a]);
//    temp.push_back(ielem_mix[0][i*9 + b]);
//    temp.push_back(ielem_mix[0][i*9 + c]);
//    temp.push_back(ielem_mix[0][i*9 + d]);


}

std::sort( temp.begin(), temp.end());
long long int hash_int;
int  t;
hash_int =0;
t = 0;
int factor;

factor = ceil(log10(n_cells));


for(std::vector<int>::iterator it = temp.begin(); it != temp.end(); ++it) {

   hash_int = hash_int + *it * pow(10,factor * t);
   t++;
}

face_labels.push_back(hash_int);


}



void unstructured_mesh::calc_face_neighbours(global_variables &globals){

ostringstream strs;

string temp;


std::fill_n( n_node, total_cells, -1 );
std::fill_n( e_node, total_cells, -1 );
std::fill_n( w_node, total_cells, -1 );
std::fill_n( s_node, total_cells, -1 );
std::fill_n( f_node, total_cells, -1 );
std::fill_n( b_node, total_cells, -1 );

// get face_labels

    for (int i = 0; i< n_cells; i++){
            //get the north and east flux streaming time step
           // west = 0, south = 1;  east = 2; north = 3; front = 4; back =5;

            //west
            face_neighbours_hash(1,5,8,4,i,globals);

            //south
            face_neighbours_hash(2,1,5,6,i,globals);

           //east
            face_neighbours_hash(2,6,7,3,i,globals);

           //north
            face_neighbours_hash(3,4,8,7,i,globals);

           //front
           face_neighbours_hash(5,6,7,8,i,globals);

           //back
           face_neighbours_hash(1,2,3,4,i,globals);



    }

    bool skip;
    int k;
    vector<long long int>::iterator it ,it1;
    int p,l;
     int element;


    std::ofstream face_labels_log;
    std::string face_labels_log_dir;
    face_labels_log_dir = globals.output_file +"/face_labels.txt";

    face_labels_log.open(face_labels_log_dir.c_str(), ios::out);

      for ( int index = 0; index < face_labels.size(); index++ ){

        face_labels_log << face_labels[index] << endl;

      }
    face_labels_log.close();

    //loop through each face
    for ( int index = 0; index < face_labels.size(); index++){

        skip = false;

        //check if found already
        k = index%6;
        element = floor(index/6);

         switch(k) {
            case 0: // west

                if(w_node[element] > -1){
                    skip = true;
                }
                 break;
            case 1: // south

                if(s_node[element] > -1){
                    skip = true;
                }
                 break;
            case 2: // east

                if(e_node[element] > -1){
                    skip = true;
                }
                 break;
            case 3: // north

                if(n_node[element] > -1){
                    skip = true;
                }
                 break;
            case 4: // front

                if(f_node[element] > -1){
                    skip = true;
                }
                 break;
            case 5: // back

                if(b_node[element] > -1){
                    skip = true;
                }
                 break;

        }

        if(skip){
            //next iteration of loop
        }else{
            //find operator
            it = std::find(face_labels.begin() +index +1,face_labels.end(),face_labels[index]);
            if( it!= face_labels.end()){
                p = std::distance( face_labels.begin(), it );

                populate_face_neighbour(index,p);
            }
        }

    }
    return;
}


void unstructured_mesh::populate_face_neighbour(int index1, int index2){

    int element1,element2,k;

    element1 = floor(index1/6);
    element2 = floor(index2/6);

    k = index1%6;



    switch(k) {
            case 0: // west
                w_node[element1] = element2;
                 break;
            case 1: // south
                s_node[element1] = element2;
                 break;
            case 2: // east
                e_node[element1] = element2;
                 break;
            case 3: // north
                n_node[element1] = element2;
                 break;
            case 4: // front
                f_node[element1] = element2;
                 break;
            case 5: // back
                b_node[element1] = element2;
                 break;
        }

    k = index2%6;

    switch(k) {
            case 0: // west
                w_node[element2] = element1;
                break;
            case 1: // south
                s_node[element2] = element1;
                break;
            case 2: // east
                e_node[element2] = element1;
                break;
            case 3: // north
                n_node[element2] = element1;
                break;
            case 4: // front
                f_node[element2] = element1;
                break;
            case 5: // back
                b_node[element2] = element1;
                break;
        }




}

void unstructured_mesh::initialise_mesh_variables(){

    total_cells = n_cells + num_bc_cells;


    centroid_x = new double [total_cells +1]();
        if (centroid_x==NULL) exit (1);
    centroid_y = new double [total_cells +1]();
        if (centroid_y==NULL) exit (1);
    centroid_z = new double [total_cells +1]();
        if (centroid_z==NULL) exit (1);

    cell_volume = new double [total_cells +1]();
        if (cell_volume==NULL) exit (1);

    total_cell_area = new double [total_cells +1]();
        if (total_cell_area==NULL) exit (1);



    cell_cfl_x = new double [total_cells +1]();
        if (cell_cfl_x==NULL) exit (1);

    cell_cfl_y = new double [total_cells +1]();
        if (cell_cfl_y==NULL) exit (1);

    cell_cfl_z = new double [total_cells +1]();
        if (cell_cfl_z==NULL) exit (1);

    cell_cfl_r = new double [total_cells +1]();
        if (cell_cfl_r==NULL) exit (1);


    for( int i =0; i< total_cells; i++){
        cell_volume[i] =0.0;
        centroid_x[i] = 0.0;
        centroid_y[i] = 0.0;
        centroid_z[i] = 0.0;
        total_cell_area[i] =0.0;
        cell_cfl_x[i] = 0.0;
        cell_cfl_y[i] = 0.0;
        cell_cfl_z[i] = 0.0;
        cell_cfl_r[i] = 1000000.0;
        }


    // face variables


    face_x = new double [n_faces +1]();
        if (face_x==NULL) exit (1);
    face_y = new double [n_faces +1]();
        if (face_y==NULL) exit (1);
    face_z = new double [n_faces +1]();
        if (face_z==NULL) exit (1);

     // area of faces
    face_area = new double [n_faces +1]();
        if (face_area==NULL) exit (1);

      delta_t_face = new double [n_faces +1]();
        if (delta_t_face==NULL) exit (1);


    for( int i =0; i< n_faces; i++){
        face_area[i] =0.0;
        face_x[i] = 0.0;
        face_y[i] = 0.0;
        face_z[i] = 0.0;
        delta_t_face[i] = 1;
    }

    //face normal
    face_i = new double [n_faces +1]();
        if (face_i==NULL) exit (1);
    face_j = new double [n_faces +1]();
        if (face_j==NULL) exit (1);
    face_k = new double [n_faces +1]();
        if (face_k==NULL) exit (1);



    //centre nodes of faces
    north_x = new double [total_cells +1]();
        if (north_x==NULL) exit (1);
    north_y = new double [total_cells +1]();
        if (north_y==NULL) exit (1);
    north_z = new double [total_cells +1]();
        if (north_z==NULL) exit (1);
    east_x = new double [total_cells +1]();
        if (east_x==NULL) exit (1);
    east_y = new double [total_cells +1]();
        if (east_y==NULL) exit (1);
    east_z = new double [total_cells +1]();
        if (east_z==NULL) exit (1);
    west_x = new double [total_cells +1]();
        if (west_x==NULL) exit (1);
    west_y = new double [total_cells +1]();
        if (west_y==NULL) exit (1);
    west_z = new double [total_cells +1]();
        if (west_z==NULL) exit (1);
    south_x = new double [total_cells +1]();
        if (south_x==NULL) exit (1);
    south_y = new double [total_cells +1]();
        if (south_y==NULL) exit (1);
    south_z = new double [total_cells +1]();
        if (south_z==NULL) exit (1);
    front_x = new double [total_cells +1]();
        if (front_x==NULL) exit (1);
    front_y = new double [total_cells +1]();
        if (front_y==NULL) exit (1);
    front_z = new double [total_cells +1]();
        if (front_z==NULL) exit (1);
    back_x = new double [total_cells +1]();
        if (back_x==NULL) exit (1);
    back_y = new double [total_cells +1]();
        if (back_y==NULL) exit (1);
    back_z = new double [total_cells +1]();
        if (back_z==NULL) exit (1);


    // area of faces
    n_area = new double [total_cells +1]();
        if (n_area==NULL) exit (1);
    s_area = new double [total_cells +1]();
        if (s_area==NULL) exit (1);
    w_area = new double [total_cells +1]();
        if (w_area==NULL) exit (1);
    e_area = new double [total_cells +1]();
        if (e_area==NULL) exit (1);
    f_area = new double [total_cells +1]();
        if (f_area==NULL) exit (1);
    b_area = new double [total_cells +1]();
        if (b_area==NULL) exit (1);




    //
    n_i = new double [total_cells +1]();
        if (n_i==NULL) exit (1);
    n_j = new double [total_cells +1]();
        if (n_j==NULL) exit (1);
    n_k = new double [total_cells +1]();
        if (n_k==NULL) exit (1);
     e_i = new double [total_cells +1]();
        if (e_i==NULL) exit (1);
    e_j = new double [total_cells +1]();
        if (e_j==NULL) exit (1);
    e_k = new double [total_cells +1]();
        if (e_k==NULL) exit (1);
     w_i = new double [total_cells +1]();
        if (w_i==NULL) exit (1);
    w_j = new double [total_cells +1]();
        if (w_j==NULL) exit (1);
    w_k = new double [total_cells +1]();
        if (w_k==NULL) exit (1);
     s_i = new double [total_cells +1]();
        if (s_i==NULL) exit (1);
    s_j = new double [total_cells +1]();
        if (s_j==NULL) exit (1);
    s_k = new double [total_cells +1]();
        if (s_k==NULL) exit (1);
    f_i = new double [total_cells +1]();
        if (s_i==NULL) exit (1);
    f_j = new double [total_cells +1]();
        if (s_j==NULL) exit (1);
    f_k = new double [total_cells +1]();
        if (s_k==NULL) exit (1);
    b_i = new double [total_cells +1]();
        if (s_i==NULL) exit (1);
    b_j = new double [total_cells +1]();
        if (s_j==NULL) exit (1);
    b_k = new double [total_cells +1]();
        if (s_k==NULL) exit (1);


    n_node = new int [total_cells +1]();
        if (n_node==NULL) exit (1);
    s_node = new int [total_cells +1]();
        if (s_node==NULL) exit (1);
    w_node = new int [total_cells +1]();
        if (w_node==NULL) exit (1);
    e_node = new int [total_cells +1]();
        if (e_node==NULL) exit (1);
    f_node = new int [total_cells +1]();
        if (f_node==NULL) exit (1);
    b_node = new int [total_cells +1]();
        if (b_node==NULL) exit (1);


    calcs_per_cell = new int [total_cells +1]();
        if (calcs_per_cell==NULL) exit (1);



    delta_t = new double [total_cells +1]();
        if (delta_t==NULL) exit (1);
    delta_t_n = new double [total_cells +1]();
    if (delta_t==NULL) exit (1);
    delta_t_e = new double [total_cells +1]();
    if (delta_t==NULL) exit (1);

     bc_neighbour = new int [num_bc_cells +1]();
    if (bc_neighbour==NULL) exit (1);

    bc_face_x = new double [num_bc_cells +1]();
        if (bc_face_x==NULL) exit (1);
    bc_face_y = new double [num_bc_cells +1]();
        if (bc_face_y==NULL) exit (1);
    bc_face_z = new double [num_bc_cells +1]();
        if (bc_face_z==NULL) exit (1);

    bc_face_direction = new double [num_bc_cells +1]();
        if (bc_face_direction==NULL) exit (1);

}



void unstructured_mesh::import_openfoam_BC_tags(global_variables &globals){

    std::string input_file;
    input_file = globals.import_file + "boundary";

    std::ifstream boundary;


    std::string line;
    std::string token1,token2;
    std:: string bc_type;
    int t,i,num_bc_tags,n_bc_faces;
    t = 0;
    i =0;
	n_wall_cells = 0;
    char comp;
    bool start = false;
    //read in points file
    boundary.open(input_file);

    while (std::getline(boundary, line))
    {

        std::istringstream iss(line);


        if(line.length() >0){
            comp = line.at(0);
        }



        if(comp == ')'){
            start = false;
        }

        if(start && i < n_faces){
            iss >> token1 >> token2;

            if (token1 == "type"){
                token2.pop_back();

                bc_type = token2;

            }

            if( token1 == "nFaces"){
                token2.pop_back();
                n_bc_faces = std::stoi(token2);

                 for( int j = 0; j< (int)n_bc_faces; j++){
                    bc_types.push_back(bc_type);
                }
				if (bc_type == "wall" || bc_type == "Wall") {
					 n_wall_cells = n_wall_cells + n_bc_faces;
				}

            }
        }

        //order of check is important
         if(comp == '('){
            start = true;
        }





        t++;
        // process pair (a,b)
    }

    boundary.close();

}


void unstructured_mesh::tecplot_output_polyhedral_unstructured_grid(global_variables &globals)
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
    double timestamp = 0.0;
    tt << timestamp;
    int fileType_e = 1;


    output_location = globals.output_file + "/plt/grid/grid.plt";
    valueLocation = new int[3];
    for(int i = 0; i < 3;i++){
           valueLocation[i] = 1;
    }
    strandID  = 0;   /* StaticZone */

   //enum fileType_e { FULL = 0, GRID = 1, SOLUTION = 2 };

     float *nodx, *nody, *nodz;
    int *connectivity;
    double solTime;
    INTEGER4 debug, i, j, k, dIsDouble, vIsDouble, zoneType,  parentZn, isBlock;
    INTEGER4 iCellMax, jCellMax, kCellMax, nFConns, fNMode, shrConn, fileType;


    INTEGER4 fileFormat; // 0 == PLT, 1 == SZPLT
    fileFormat = 0;

    INTEGER4 nNodes, nCells, nFaces, connectivityCount, index;

    int t,r;
    int n_ghost;  // number of host cells;

    n_ghost = 1;

    if(globals.testcase == 3){
        n_ghost = 2;
    }

    debug     = 1;
    vIsDouble = 0;
    dIsDouble = 0;
    nNodes = n_vertices;
    nCells =  n_cells;
    nFaces = n_faces; /* Not used */
    zoneType  = 7;      /* FEPolyhedron */
    solTime   = timestamp;

    parentZn  = 0;      /* No Parent */
    isBlock   = 1;      /* Block */
    iCellMax  = 0;
    jCellMax  = 0;
    kCellMax  = 0;
    nFConns   = 0;
    fNMode    = 1;
    shrConn   = 0;
    fileType  = fileType_e;
    ss << nCells;
    zone_name = ss.str();

    /* The number of face nodes in the zone. This example creates
    * a zone with a single pyramidal cell. This cell has four
    * triangular faces and one rectangular face, yielding a total
    * of 16 face nodes.
    */
    INTEGER4 NumFaceNodes = n_face_points;
    INTEGER4 NumBConns  = 0;    /* No Boundary Connections */
    INTEGER4 NumBItems = 0;     /* No Boundary Items */


    /*
     * Open the file and write the tecplot datafile
     * header information
     */
          i = TECINI142((char*) "Polyhedral Mesh" ,
                  (char*)"nodx nody nodz",
                  (char*) output_location.c_str(),
                  (char*) ".",
                  &fileFormat,
                  &fileType,
                  &debug,
                  &vIsDouble);




    i = TECAUXSTR142("Re" ,  reynolds_text.c_str());


    nodx  = (float*)calloc(nNodes , sizeof(float));
    nody  = (float*)calloc(nNodes , sizeof(float));
    nodz  = (float*)calloc(nNodes , sizeof(float));
    t = 0;
    r=0;



    std::string filename;
    std::ofstream globals_txt,globals_txt1;
    std::string globals_file,globals_file1;
    output_location = globals.output_file;
    filename = globals.simulation_name;
    globals_file = output_location + "/trial.txt";
    globals_file1 = output_location + "/trial1.txt";
    globals_txt.open(globals_file.c_str(), ios::out);


    for (k = 0; k < n_vertices; k++){

          nodx[t] = (float)(x[t]);
          nody[t] = (float)(y[t]);
          nodz[t] = (float)(z[t]);

          globals_txt << t << "," << nodx[t] << ","<< nody[t] << ","<< nodz[t] << std::endl;
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


      INTEGER4 *FaceNodeCounts = new INTEGER4[n_faces];
      INTEGER4 *FaceNodes = new INTEGER4[NumFaceNodes];
      int q;
      q =0;
      for(int k =0; k < n_faces; k++){
        FaceNodeCounts[k] =face_connectivity[k].size();
            for(int l =0;l<face_connectivity[k].size();l++){
                FaceNodes[q] = face_connectivity[k][l] +1;
                q++;
            }
      }

      INTEGER4 *FaceLeftElems = new INTEGER4[n_faces];
      INTEGER4 *FaceRightElems = new INTEGER4[n_faces];
      for(k =0; k < n_faces; k++){
        FaceLeftElems[k] = owner [k] +1;
        if(neighbour[k] == 0){
            FaceRightElems[k] = 0;
        }else{
            FaceRightElems[k] = neighbour[k] +1;
        }
      }



    i = tecpolyface142(&n_faces, FaceNodeCounts,FaceNodes, FaceLeftElems,FaceRightElems);

    free(nodx);
    free(nody);
    free(nodz);

    delete [] FaceLeftElems;
    FaceLeftElems = NULL;
    delete []FaceRightElems;
    FaceRightElems = NULL;
    delete [] FaceNodes;
    FaceNodes = NULL;
    delete [] FaceNodeCounts;
    FaceNodeCounts = NULL;


    delete [] valueLocation;
    valueLocation= NULL;

    i = TECEND142();
}



void unstructured_mesh::tecplot_output_unstructured_grid(global_variables &globals)
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
    double timestamp = 0.0;
    tt << timestamp;
    int fileType_e = 1;


    output_location = globals.output_file + "/plt/grid/grid.szplt";
    valueLocation = new int[3];
    for(int i = 0; i < 3;i++){
           valueLocation[i] = 1;
    }
    strandID  = 0;   /* StaticZone */

   //enum fileType_e { FULL = 0, GRID = 1, SOLUTION = 2 };

     float *nodx, *nody, *nodz;
    int *connectivity;
    double solTime;
    INTEGER4 debug, i, j, k, dIsDouble, vIsDouble, zoneType,  parentZn, isBlock;
    INTEGER4 iCellMax, jCellMax, kCellMax, nFConns, fNMode, shrConn, fileType;


    INTEGER4 fileFormat; // 0 == PLT, 1 == SZPLT
    fileFormat = 1;

    INTEGER4 nNodes, nCells, nFaces, connectivityCount, index;

    int t,r;
    int n_ghost;  // number of host cells;

    n_ghost = 1;

    if(globals.testcase == 3){
        n_ghost = 2;
    }

    debug     = 1;
    vIsDouble = 0;
    dIsDouble = 0;
    nNodes = n_vertices;
    nCells =  n_cells;
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
          i = TECINI142((char*) "Couette Flow" ,
                  (char*)"nodx nody nodz",
                  (char*) output_location.c_str(),
                  (char*) ".",
                  &fileFormat,
                  &fileType,
                  &debug,
                  &vIsDouble);




    i = TECAUXSTR142("Re" ,  reynolds_text.c_str());


    nodx  = (float*)calloc(nNodes , sizeof(float));
    nody  = (float*)calloc(nNodes , sizeof(float));
    nodz  = (float*)calloc(nNodes , sizeof(float));
    t = 0;
    r=0;



    std::string filename;
    std::ofstream globals_txt,globals_txt1;
    std::string globals_file,globals_file1;
    output_location = globals.output_file;
    filename = globals.simulation_name;
    globals_file = output_location + "/trial.txt";
    globals_file1 = output_location + "/trial1.txt";
    globals_txt.open(globals_file.c_str(), ios::out);


    for (k = 0; k < n_vertices; k++){

          nodx[t] = (float)(x[t])/X;
          nody[t] = (float)(y[t])/Y;
          nodz[t] = (float)(z[t])/Z;


          globals_txt << t << "," << nodx[t] << ","<< nody[t] << ","<< nodz[t] << std::endl;
          t++;

        }


      globals_txt.close();
      globals_txt1.open(globals_file1.c_str(), ios::out);
    connectivityCount = 8 * nCells;
    connectivity = (INTEGER4*)malloc(connectivityCount * sizeof(INTEGER4));
    t =0;

    if (globals.import_format == "openfoam"){
        for (int k = 0; k < (n_cells ); k++)

            {
                index = k  * 8;
                t = k*9;
                connectivity[index] = volume_connectivity[k][1] +1;
                connectivity[index + 1] = volume_connectivity[k][5] +1;
                connectivity[index + 2] = volume_connectivity[k][4]+1;
                connectivity[index + 3] = volume_connectivity[k][0]+1;
                connectivity[index + 4] = volume_connectivity[k][3]+1;
                connectivity[index + 5] = volume_connectivity[k][7]+1;
                connectivity[index + 6] = volume_connectivity[k][6]+1;
                connectivity[index + 7] = volume_connectivity[k][2]+1;


                globals_txt1 << connectivity[index] << ","<< connectivity[index + 1] << ","<< connectivity[index + 2] << ","<< connectivity[index + 3] << ","
                << connectivity[index + 4] << ","<< connectivity[index + 5] << ","<< connectivity[index + 6] << ","<<connectivity[index + 7] << std::endl;

           }

    }else{
        for (int k = 0; k < (n_cells ); k++)

            {
                index = k  * 8;
                t = k*9;
//                connectivity[index] = ielem_mix[0][t+1];
//                connectivity[index + 1] = ielem_mix[0][t+2];
//                connectivity[index + 2] = ielem_mix[0][t+3];
//                connectivity[index + 3] = ielem_mix[0][t+4];
//                connectivity[index + 4] = ielem_mix[0][t+5];
//                connectivity[index + 5] = ielem_mix[0][t+6];
//                connectivity[index + 6] = ielem_mix[0][t+7];
//                connectivity[index + 7] = ielem_mix[0][t+8];
//
//
//                globals_txt1 << connectivity[index] << ","<< connectivity[index + 1] << ","<< connectivity[index + 2] << ","<< connectivity[index + 3] << ","
//                << connectivity[index + 4] << ","<< connectivity[index + 5] << ","<< connectivity[index + 6] << ","<<connectivity[index + 7] << std::endl;

                   }


        globals_txt1.close();


    }
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

    i = TECDAT142(&nNodes, nodx, &dIsDouble);
    i = TECDAT142(&nNodes, nody, &dIsDouble);
    i = TECDAT142(&nNodes, nodz, &dIsDouble);
   i = TECNODE142(&connectivityCount, connectivity);
    free(connectivity);
    free(nodx);
    free(nody);
    free(nodz);


    delete [] valueLocation;
    valueLocation= NULL;

    i = TECEND142();
}




void unstructured_mesh::import_tecplot_mesh(global_variables &globals){

    std::string input_file;
    input_file = globals.import_file + ".szplt";


    // Open a .szplt file for reading with TECIO
    void* inputFileHandle = NULL;
    int i = tecFileReaderOpen(input_file.c_str(), &inputFileHandle);

    // Read the characteristics of the data set
    char* dataSetTitle = NULL;
    i = tecDataSetGetTitle(inputFileHandle, &dataSetTitle);

    int32_t numVars;
    i = tecDataSetGetNumVars(inputFileHandle, &numVars);
    std::ostringstream outputStream;
    for (int32_t var = 1; var <= numVars; ++var)
    {
        char* name = NULL;
        i = tecVarGetName(inputFileHandle, var, &name);
        outputStream << name;
        if (var < numVars)
            outputStream << ',';
        tecStringFree(&name);
    }

    int32_t fileType;
    i = tecFileGetType(inputFileHandle, &fileType);

    int32_t numZones;
    i = tecDataSetGetNumZones(inputFileHandle, &numZones);

    int32_t isDouble = 0;
    int32_t const FieldDataType_Double = 2; // ...TecUtil types are not available to TecIO so redefine
    if (numZones > 0)
    {
        int32_t type;
        i = tecZoneVarGetType(inputFileHandle, 1, 1, &type);
        if (type == FieldDataType_Double)
            isDouble = 1;
    }

  for (int32_t inputZone = 1; inputZone <= numZones; ++inputZone)
    {
        // Retrieve zone characteristics
        int32_t zoneType;
        i = tecZoneGetType(inputFileHandle, inputZone, &zoneType);
        if (zoneType == 6 || zoneType == 7)
            throw std::runtime_error("Unsupported zone type.");

        char* zoneTitle = NULL;
        i = tecZoneGetTitle(inputFileHandle, inputZone, &zoneTitle);

        int64_t iMax, jMax, kMax;
        i = tecZoneGetIJK(inputFileHandle, inputZone, &iMax, &jMax, &kMax);

        std::vector<int32_t> varTypes(numVars);
        std::vector<int32_t> passiveVarList(numVars);
        std::vector<int32_t> valueLocation(numVars);
        std::vector<int32_t> shareVarFromZone(numVars);
        for (int32_t var = 1; var <= numVars; ++var)
        {
            i = tecZoneVarGetType(inputFileHandle, inputZone, var, &varTypes[var - 1]);
            i = tecZoneVarGetSharedZone(inputFileHandle, inputZone, var, &shareVarFromZone[var - 1]);
            i = tecZoneVarGetValueLocation(inputFileHandle, inputZone, var, &valueLocation[var - 1]);
            i = tecZoneVarIsPassive(inputFileHandle, inputZone, var, &passiveVarList[var - 1]);
        }

        int32_t shareConnectivityFromZone;
        i = tecZoneConnectivityGetSharedZone(inputFileHandle, inputZone, &shareConnectivityFromZone);

        int32_t faceNeighborMode;
        i = tecZoneFaceNbrGetMode(inputFileHandle, inputZone, &faceNeighborMode);

        int64_t numFaceConnections;
        i = tecZoneFaceNbrGetNumConnections(inputFileHandle, inputZone, &numFaceConnections);

        int64_t numFaceConnections1;
        i = tecZoneFaceNbrGetNumValues(inputFileHandle, inputZone, &numFaceConnections1);


        double solutionTime;
        int32_t strandID;
        i = tecZoneGetSolutionTime(inputFileHandle, inputZone, &solutionTime);
        i = tecZoneGetStrandID(inputFileHandle, inputZone, &strandID);

        int32_t parentZone;
        i = tecZoneGetParentZone(inputFileHandle, inputZone, &parentZone);

        tecStringFree(&zoneTitle);

        // Retrieve zone data and send to the output file
        for (int32_t var = 1; var <= numVars; ++var)
        {
            if (passiveVarList[var - 1] == 0 && shareVarFromZone[var - 1] == 0)
            {
                int64_t numValues;
                i = tecZoneVarGetNumValues(inputFileHandle, inputZone, var, &numValues);
                // For large zones, could "chunk" this input/output--read/write the var in pieces instead of all at once
                switch((FieldDataType_e)varTypes[var - 1])
                {
                case FieldDataType_Float:
                    {
                        std::vector<float> values(numValues);
                        i = tecZoneVarGetFloatValues(inputFileHandle, inputZone, var, 1, numValues, &values[0]);
                        //i = tecZoneVarWriteFloatValues(outputFileHandle, outputZone, var, 0, numValues, &values[0]);
                    }
                    break;
                case FieldDataType_Double:
                    {
                        std::vector<double> values(numValues);
                        i = tecZoneVarGetDoubleValues(inputFileHandle, inputZone, var, 1, numValues, &values[0]);
                       // i = tecZoneVarWriteDoubleValues(outputFileHandle, outputZone, var, 0, numValues, &values[0]);
                    }
                    break;
                case FieldDataType_Int32:
                    {
                        std::vector<int32_t> values(numValues);
                        i = tecZoneVarGetInt32Values(inputFileHandle, inputZone, var, 1, numValues, &values[0]);
                        //i = tecZoneVarWriteInt32Values(outputFileHandle, outputZone, var, 0, numValues, &values[0]);
                    }
                    break;
                case FieldDataType_Int16:
                    {
                        std::vector<int16_t> values(numValues);
                        i = tecZoneVarGetInt16Values(inputFileHandle, inputZone, var, 1, numValues, &values[0]);
                        //i = tecZoneVarWriteInt16Values(outputFileHandle, outputZone, var, 0, numValues, &values[0]);
                    }
                    break;
                case FieldDataType_Byte:
                    {
                        std::vector<uint8_t> values(numValues);
                        i = tecZoneVarGetUInt8Values(inputFileHandle, inputZone, var, 1, numValues, &values[0]);
                        //i = tecZoneVarWriteUInt8Values(outputFileHandle, outputZone, var, 0, numValues, &values[0]);
                    }
                    break;
                default:
                    i = -1;
                    break;
                }
            }
        }

        // Write zone face neighbors, if any
        if (numFaceConnections > 0)
        {
            int64_t numFaceValues;
            i = tecZoneFaceNbrGetNumValues(inputFileHandle, inputZone, &numFaceValues);
            int32_t are64Bit;
            i = tecZoneFaceNbrsAre64Bit(inputFileHandle, inputZone, &are64Bit);
            if (are64Bit)
            {
                std::vector<int64_t> faceConnections(numFaceValues);
                i = tecZoneFaceNbrGetConnections64(inputFileHandle, inputZone, &faceConnections[0]);
                //i = tecZoneFaceNbrWriteConnections64(outputFileHandle, outputZone, &faceConnections[0]);
            }
            else
            {
                std::vector<int32_t> faceConnections(numFaceValues);
                i = tecZoneFaceNbrGetConnections(inputFileHandle, inputZone, &faceConnections[0]);
               // i = tecZoneFaceNbrWriteConnections32(outputFileHandle, outputZone, &faceConnections[0]);
            }
        }

        int64_t numValues;
        i = tecZoneNodeMapGetNumValues(inputFileHandle, inputZone, jMax, &numValues);
        std::vector<int64_t> nodeMap_64(numValues);
        std::vector<int32_t> nodeMap_32(numValues);
        // Retrieve zone node map, if any, and send to the output file
        if (zoneType != 0 && shareConnectivityFromZone == 0)
        {

            int32_t is64Bit;
            i = tecZoneNodeMapIs64Bit(inputFileHandle, inputZone, &is64Bit);
            if (is64Bit)
            {

                i = tecZoneNodeMapGet64(inputFileHandle, inputZone, 1, jMax, &nodeMap_64[0]);
              //  i = tecZoneNodeMapWrite64(outputFileHandle, outputZone, 0, 1, numValues, &nodeMap[0]);
            }
            else
            {

                i = tecZoneNodeMapGet(inputFileHandle, inputZone, 1, jMax, &nodeMap_32[0]);
               // i = tecZoneNodeMapWrite32(outputFileHandle, outputZone, 0, 1, numValues, &nodeMap[0]);
            }
        }

        // Retrieve and write any zone aux data
        int32_t numItems;
        i = tecZoneAuxDataGetNumItems(inputFileHandle, inputZone, &numItems);
        for (int32_t whichItem = 1; whichItem <= numItems; ++whichItem)
        {
            char* name = NULL;
            char* value = NULL;
            i = tecZoneAuxDataGetItem(inputFileHandle, inputZone, whichItem, &name, &value);
           // i = tecZoneAddAuxData(outputFileHandle, outputZone, name, value);
            tecStringFree(&name);
            tecStringFree(&value);
        }
    }

    int j;
    j = 1;


}
