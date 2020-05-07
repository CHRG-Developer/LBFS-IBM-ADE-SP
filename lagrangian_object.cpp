#include "lagrangian_object.h"
#include "global_variables.h"
#include <cmath>
#include <iostream>
#include <fstream> 
#include <sstream>
#include <vector>
#include <algorithm> 
#include <string>
#include "vector_var.h"



lagrangian_object::lagrangian_object()
{
}


lagrangian_object::~lagrangian_object()
{
	/*delete[] node_x_ref;
	node_x_ref = NULL;
	delete[] node_y_ref;
	node_y_ref = NULL;
	delete[] node_z_ref;
	node_z_ref = NULL;*/


}

void lagrangian_object::initialise(double PI) {

	// rigid 3D Cylinder
	if (type == 1) {
		node_x_ref = new double[num_nodes];
		if (node_x_ref == NULL) exit(1);
		node_y_ref = new double[num_nodes];
		if (node_y_ref == NULL) exit(1);
		node_z_ref = new double[num_nodes];
		if (node_z_ref == NULL) exit(1);

		node_x = new double[num_nodes];
		if (node_x == NULL) exit(1);
		node_y = new double[num_nodes];
		if (node_y == NULL) exit(1);
		node_z = new double[num_nodes];
		if (node_z == NULL) exit(1);

		populate_node_reference_displacement(PI);
	}/// spring network mesh
	else if (type == 2) {
		import_network_mesh();

	}


	node_vel_x = new double[num_nodes];
	if (node_vel_x == NULL) exit(1);
	node_vel_y = new double[num_nodes];
	if (node_vel_y == NULL) exit(1);
	node_vel_z = new double[num_nodes];
	if (node_vel_z == NULL) exit(1);

	node_force_x = new double[num_nodes];
	if (node_force_x == NULL) exit(1);
	node_force_y = new double[num_nodes];
	if (node_force_y == NULL) exit(1);
	node_force_z = new double[num_nodes];
	if (node_force_z == NULL) exit(1);

	external_loading_factor = new double[num_nodes];
	if (external_loading_factor == NULL) exit(1);

	//num_springs = 0;
	//Mesh_Process();
	compatibility();

	min_index = 0;
	max_index = 0;


	spontaneous_angle = sqrt(3) * (num_nodes - 2) - 5 * PI;
	spontaneous_angle = spontaneous_angle / (sqrt(3) * (num_nodes - 2) - 3 * PI);
	spontaneous_angle = acos(spontaneous_angle);

}

void lagrangian_object::apply_optical_tweezers_force() {


	// get 2% of nodes
	 int two_percent = ceil(num_nodes * 0.02);
	//two_percent = 1;
	 boundary_force = boundary_force / two_percent;

	// get min x nodes - 2% of nodes

	std::vector<int> index;
	std::vector<double> magnitude;
	

	for (int i = 0; i < num_nodes; i++) {

		//use loop to initialise loading while we're at it
		external_loading_factor[i] = 0.0;

		//fill vectors
		if (i < two_percent ) {
			index.push_back(i);
			magnitude.push_back(sorted_coordinate[i*3]);
		}
		else {
			//check if smaller than min
			max_index = std::max_element(magnitude.begin(), magnitude.end()) - magnitude.begin();

			if (sorted_coordinate[i * 3] < magnitude[max_index]) {
				magnitude[max_index] = sorted_coordinate[i * 3];
				index[max_index] = i;
			}
		}
	}
	
	for (int i = 0; i < index.size(); i++) {
		external_loading_factor[index[i]] = -1.0;
	}



	index.clear();
	magnitude.clear();


	// get max x node - 2% of nodes

	for (int i = 0; i < num_nodes; i++) {

	
		//fill vectors
		if (i < two_percent ) {
			index.push_back(i);
			magnitude.push_back(sorted_coordinate[i * 3]);
		}
		else {
			//check if smaller than min
			min_index = std::min_element(magnitude.begin(), magnitude.end()) - magnitude.begin();

			if (sorted_coordinate[i * 3] > magnitude[min_index]) {
				magnitude[min_index] = sorted_coordinate[i * 3];
				index[min_index] = i;
			}
		}
	}

	for (int i = 0; i < index.size(); i++) {
		external_loading_factor[index[i]] = 1.0;
	}
		
	 //finally store index of min and max index for later calculations

	double max_val, min_val;
	max_val = -1000000000000;
	min_val = 100000000000;
	for (int i = 0; i < num_nodes; i++) {
		if (sorted_coordinate[i * 3] > max_val) {
			max_index = i;
			max_val = sorted_coordinate[i * 3];

		}
		if (sorted_coordinate[i * 3] < min_val) {
			min_index = i;
			min_val = sorted_coordinate[i * 3];

		}
	}

}

void lagrangian_object::import_network_mesh() {
	std::cout << "Constructing a mesh... " << std::endl;

	//mesh_read(mesh_file);
	// read the stress free mesh originally
	// connectivity etc. is the same, update the xyz nodes after
	//

	if (mz_importer) {
		mesh_read(stress_free_mesh);
	}
	else {

		osay_mesh_read(stress_free_mesh);

	}


	/*get_nod_max_freedom();
	get_nod_min_freedom();
*/
	std::cout << "\tSuccessfully constructed the mesh!\n" << std::endl;

}
void lagrangian_object::populate_node_reference_displacement(double PI) {

	//rigid cylinder in x direction fully 3d
	if (type == 1) {


		depth_nodes = (int)sqrt(num_nodes);
		radial_nodes = depth_nodes;

		double delta_z = depth / depth_nodes;
		int k = 0;


		for (int j = 0; j < depth_nodes; j++) {
			for (int i = 0; i < radial_nodes; i++) {

				node_x_ref[k] = centre_x + radius * sin(2 * PI * i / radial_nodes);
				node_y_ref[k] = centre_y + radius * cos(2 * PI * i / radial_nodes);
				node_z_ref[k] = centre_z + depth / 2 - j * depth / depth_nodes - delta_z / 2;

				node_x[k] = centre_x + radius * sin(2 * PI * i / radial_nodes);
				node_y[k] = centre_y + radius * cos(2 * PI * i / radial_nodes);
				node_z[k] = centre_z + depth / 2 - j * depth / depth_nodes - delta_z / 2;;

				k = k + 1;

			}
		}

		delta_z = 1;

	}


}

void lagrangian_object::osay_mesh_read(std::string inputFileName) {

	std::cout << "\tReading mesh file..." << std::endl;

	//read in the mesh file.
	std::ifstream inputFile;
	inputFile.open(inputFileName.c_str());
	if (inputFile.is_open()) {

		//std::cout << "Successfully opened the NODES input file! " << std::endl;

		//for nodes.
		inputFile >> num_nodes;
		inputFile >> num_tets;

		//std::cout << num_nodes << " " << triN << std::endl;

		//update x_y_z nodes
		node_x = new double[num_nodes];
		if (node_x == NULL) exit(1);
		node_y = new double[num_nodes];
		if (node_y == NULL) exit(1);
		node_z = new double[num_nodes];
		if (node_z == NULL) exit(1);
		node_x_ref = new double[num_nodes];
		if (node_x_ref == NULL) exit(1);
		node_y_ref = new double[num_nodes];
		if (node_y_ref == NULL) exit(1);
		node_z_ref = new double[num_nodes];
		if (node_z_ref == NULL) exit(1);

		//import xyz nodes
		for (int i = 0; i < num_nodes; i++) {

			inputFile >> node_x[i];
			inputFile >> node_y[i];
			inputFile >> node_z[i];
			node_x_ref[i] = node_x[i] + centre_x;
			node_y_ref[i] = node_y[i] + centre_y;
			node_z_ref[i] = node_z[i] + centre_z;
			node_x[i] = node_x_ref[i];
			node_y[i] = node_y_ref[i];
			node_z[i] = node_z_ref[i];

		}
		//get num_nodes
		num_node_neighbours = new int[num_nodes];

		//use vectors to collate Node connectivity as it's unstructured
		for (int i = 0; i < num_nodes; i++) {
			node_neighbours.push_back(std::vector <int>());
			num_node_neighbours[i] = 0;
		}
		//get tet connectivity and degree of freedom of nodes
		int tet = 0;
		tet_connectivity = new int[num_tets * 3]; {
			for (int i = 0; i < num_tets * 3; i++) {
				int k;
				inputFile >> k;
				tet_connectivity[i] = k - 1;

				num_node_neighbours[k - 1]++; //add connectivity for each tet connectivity. One spring for each triangle

			}
		}

		int node_centre, pair1, pair2;
		//push back every tet pair
		for (int i = 0; i < num_tets ; i++) {
			
			node_centre = tet_connectivity[i*3 + 0];
			pair1 = tet_connectivity[i * 3 + 1];
			pair2 = tet_connectivity[i * 3 + 2];

			node_neighbours[node_centre].push_back(pair1);
			node_neighbours[node_centre].push_back(pair2);

			node_centre = tet_connectivity[i * 3 + 1];
			pair1 = tet_connectivity[i * 3 + 2];
			pair2 = tet_connectivity[i * 3 + 0 ];

			node_neighbours[node_centre].push_back(pair1);
			node_neighbours[node_centre].push_back(pair2);

			node_centre = tet_connectivity[i * 3 + 2];
			pair1 = tet_connectivity[i * 3 + 0];
			pair2 = tet_connectivity[i * 3 + 1];

			node_neighbours[node_centre].push_back(pair1);
			node_neighbours[node_centre].push_back(pair2);

		}

		//for each node -> need to sort into outward normals and sequential ordr
		vector_var p1, p2, node_c, centre;

		centre.set_equal(centre_x, centre_y, centre_z);
		int k;

		for (int i = 0; i < num_nodes; i++) {

			std::cout << i << "," << endl;

			
			//only one check of outward normal needed
			pair1 = node_neighbours[i][0];
			pair2 = node_neighbours[i][1];

			p1.set_equal(node_x[pair1], node_y[pair1], node_z[pair1]);
			p2.set_equal(node_x[pair2], node_y[pair2], node_z[pair2]);
			node_c.set_equal(node_x[i], node_y[i], node_z[i]);

			p1.subtract(node_c);
			p2.subtract(node_c);

			p1.cross_product(p2);

			//p1 is outward normal, now dot product with outward vector from centroid
			// does not work with biconcave disk
			node_c.subtract(centre);
			p1.Dot_Product(node_c);
			
			if (p1.Magnitude() < 0.0) {
				std::iter_swap(node_neighbours[i].begin(), node_neighbours[i].begin() + 1);
			}

			for (auto it = node_neighbours[i].begin(); it != node_neighbours[i].end(); it++)
			{
				std::cout << *it << ' ' ;
			}
			std::cout<<  endl;

			// first triangle is now outward normal, find next pair
			auto it = node_neighbours[i].begin() +1;

			k = 0;
			//loop untill end of connectivity
			for ( int t = 3; t < node_neighbours[i].size(); t+=2){
				//find next triangle
				auto it1 = std::find(it+1 , node_neighbours[i].end(), *it);

				int num = std::distance(it, it1);

				//swap coordinates to get in order
				if (num % 2 == 0) {
					std::iter_swap(it +1 , it +num);
					std::iter_swap(it + 2, it + num -1);
				}
				else {
					std::iter_swap(it + 1, it + num);
					std::iter_swap(it + 2, it + num + 1);
				}

				std::cout << k << "," << endl;
				//iterate to next triangle
				it = it + 2;
				k = k + 1;
			}

			for (auto it = node_neighbours[i].begin(); it != node_neighbours[i].end(); it++)
			{
				std::cout << *it << ' ';
			}
			std::cout << endl;


			//erase all duplicate values and retain order
			node_neighbours[i].erase(uniquify(node_neighbours[i].begin(), node_neighbours[i].end()), node_neighbours[i].end() );

			for (auto it = node_neighbours[i].begin(); it != node_neighbours[i].end(); it++)
			{
				std::cout << *it << ' ';
			}
			std::cout << endl;


		}
		num_springs = 0;
		for (int i = 0; i < num_nodes; i++) {
			num_springs += num_node_neighbours[i];
		}

		num_springs = num_springs / 2;

		//spring connectivity

		// will need way to remove duplicates
		spring_connectivity = new int[num_springs * 4]; 
		int a, b, c;
		std::vector <int> temp;
		std::vector <long long int> spring_hash;
		long long int hash_int;
		int t;
		int factor;
		factor = ceil(log10(num_nodes));
		int spring_index;
		spring_index = 0;

		for (int i = 0; i < num_nodes; i++) {
			std::cout << i << endl;
			for (int j =0; j < node_neighbours[i].size(); j++)
			{
				//get four nodes
				a = node_neighbours[i][j];
				b = node_neighbours[i][(j+1) % node_neighbours[i].size()];
				c = node_neighbours[i][(j + 2) % node_neighbours[i].size()];

				//get spring hash
				temp.clear();
				temp.push_back(i);
				temp.push_back(a);
				temp.push_back(b);
				temp.push_back(c);

				std::sort(temp.begin(), temp.end());

				hash_int = 0;
				t = 0;

				for (std::vector<int>::iterator it = temp.begin(); it != temp.end(); ++it) {

					hash_int = hash_int + *it * pow(10, factor * t);
					t++;
				}

				//check if it exists
				std::vector<long long int>::iterator it2 = std::find(spring_hash.begin(), spring_hash.end(), hash_int);
				
				// if it doesn't add to spring connectivity
				if (it2 == spring_hash.end()) {
					spring_hash.push_back(hash_int);
					spring_connectivity[spring_index * 4] = i;
					spring_connectivity[spring_index * 4 + 1] = a;
					spring_connectivity[spring_index * 4 + 2] = b;
					spring_connectivity[spring_index * 4 + 3] = c;
					spring_index++;
				}
				
				t = 0; // just for debuggin purposes
			}

		}

	}

	else {
		std::cout << "Failed to open the mesh input file. Press ENTER to exit!" << std::endl;
		getchar();
		exit(EXIT_FAILURE);
	}
	inputFile.close();



	return;
}

void lagrangian_object::mesh_read(std::string inputFileName) {

	std::cout << "\tReading mesh file..." << std::endl;

	//read in the mesh file.
	std::ifstream inputFile;
	inputFile.open(inputFileName.c_str());
	if (inputFile.is_open()) {

		//std::cout << "Successfully opened the NODES input file! " << std::endl;

		//for nodes.
		inputFile >> num_nodes;
		inputFile >> num_springs;
		inputFile >> num_tets;

		//std::cout << num_nodes << " " << triN << std::endl;

		//update x_y_z nodes
		node_x = new double[num_nodes];
		if (node_x == NULL) exit(1);
		node_y = new double[num_nodes];
		if (node_y == NULL) exit(1);
		node_z = new double[num_nodes];
		if (node_z == NULL) exit(1);
		node_x_ref = new double[num_nodes];
		if (node_x_ref == NULL) exit(1);
		node_y_ref = new double[num_nodes];
		if (node_y_ref == NULL) exit(1);
		node_z_ref = new double[num_nodes];
		if (node_z_ref == NULL) exit(1);

		//import xyz nodes
		for (int i = 0; i < num_nodes; i++) {

			inputFile >> node_x[i];
			inputFile >> node_y[i];
			inputFile >> node_z[i];
			node_x_ref[i] = node_x[i] + centre_x;
			node_y_ref[i] = node_y[i] + centre_y;
			node_z_ref[i] = node_z[i] + centre_z;
			node_x[i] = node_x_ref[i];
			node_y[i] = node_y_ref[i];
			node_z[i] = node_z_ref[i];

			
		}

		//get tet connectivity

		tet_connectivity = new int[num_tets * 3]; {
			for (int i = 0; i < num_tets * 3; i++) {
				int k;
				inputFile >> k;
				tet_connectivity[i] = k - 1;

			}
		}

		//spring connectivity
		spring_connectivity = new int[num_springs * 4]; {
			for (int i = 0; i < num_springs * 4; i++) {
				int k;
				inputFile >> k;
				spring_connectivity[i] = k - 1;
			}
		}

		num_node_neighbours = new int[num_nodes];
		
		//use vectors to collate Node connectivity as it's unstructured
		for (int i = 0; i < num_nodes; i++) {
			node_neighbours.push_back(std::vector <int>());
		}

		int fdm;
		int nIdx;

		for (int i = 0; i < num_nodes; i++) {

			
			inputFile >> fdm;
			num_node_neighbours[i] = fdm;

			//nod_nodIdx[i] = new int[fdm + 1];
			for (int j = 0; j < fdm; j++) {

				
				inputFile >> nIdx;
				node_neighbours[i].push_back(nIdx - 1);
			}
			
		}
	}

	else {
		std::cout << "Failed to open the mesh input file. Press ENTER to exit!" << std::endl;
		getchar();
		exit(EXIT_FAILURE);
	}
	inputFile.close();



	return;
}
void lagrangian_object::add_volume(double* vol, bool first_timestep,int n_blocks) {

	total_volume = 0.0;

	for (int i = 0; i < n_blocks; i++) {
		total_volume = total_volume + vol[i];
	}

	if (first_timestep) {
		//now account for target volume change   // 0.95 for sphere to stress free -> 0.95/0.642 for Biconcave disk
		total_volume_0 = total_volume * volume_ratio;
		//total_volume_0 = total_volume;
	}


	return;

}

void lagrangian_object::add_area(double* area, bool first_timestep, int n_blocks) {

	total_area = 0.0;

	for (int i = 0; i < n_blocks; i++) {
		total_area = total_area + area[i];
	}

	if (first_timestep) {
		total_area_0 = total_area;
	}


}


void lagrangian_object::add_MAD(double* MAD, bool first_timestep, int n_blocks) {

	total_MAD= 0.0;

	for (int i = 0; i < n_blocks; i++) {
		total_MAD = total_MAD + MAD[i];
	}
	total_MAD = total_MAD / total_area;

	if (first_timestep) {
		total_MAD_0 = total_MAD;
	}


}


void lagrangian_object::pre_process_spring_network() {

	//first step: get node neighbours

	for (int i = 0; i < num_tets; i++) {
		node_neighbours.push_back(std::vector <int>());
	}

	int k, l, m;

	for (int i = 0; i < num_tets; i++) {

		k = tet_connectivity[i * 3];
		l = tet_connectivity[i * 3 + 1];
		m = tet_connectivity[i * 3 + 2];

		node_neighbours[k].push_back(l);
		node_neighbours[k].push_back(m);
		node_neighbours[l].push_back(k);
		node_neighbours[l].push_back(m);
		node_neighbours[m].push_back(k);
		node_neighbours[m].push_back(l);

	}

	for (int i = 0; i < num_nodes; i++) {
		sort(node_neighbours[i].begin(), node_neighbours[i].end());
		node_neighbours[i].erase(unique(node_neighbours[i].begin(), node_neighbours[i].end()), node_neighbours[i].end());

	}

	int i;
	i = i;


}


/// load of shite code to go between two different importers developed by the lads
void lagrangian_object::compatibility() {

	//transfer variables from new code to old code
	// ensure no breakages and I'm lazy!

	
	nodEdge = new bool[num_nodes];
	if (nodEdge == NULL) exit(1);
	sorted_nodFreedom = new int[num_nodes];
	if (sorted_nodFreedom == NULL) exit(1);
	sorted_coordinate = new double[num_nodes * 3];

	int min, max;
	min = num_node_neighbours[0];
	max = num_node_neighbours[0];

	for (int i = 0; i < num_nodes; i++) {
		
		const int i_dfm = num_node_neighbours[i];

		if (i_dfm < min) { min = i_dfm; }

		if (i_dfm > max) { max = i_dfm; }


	}
	maxFreedom = max;
	minFreedom = min;

	nodIdx = new int[num_nodes* maxFreedom];
	if (nodIdx == NULL) exit(1);

	for (int i = 0; i < num_nodes; i++) {

		const int i_dfm = num_node_neighbours[i];


		for (int j = 0; j < i_dfm; j++) {

			nodIdx[i*maxFreedom + j] = node_neighbours[i][j];

		}
		sorted_nodFreedom[i] = i_dfm;

		sorted_coordinate[i * 3] = node_x[i];
		sorted_coordinate[i * 3 + 1] = node_y[i];
		sorted_coordinate[i * 3 + 2] = node_z[i];

	}
	
	//get tet connectivity

	sorted_triIdx = new int[num_tets * 3]; {
		for (int i = 0; i < num_tets * 3; i++) {
			sorted_triIdx[i] = tet_connectivity[i];
		}
	}

	//spring connectivity
		sprIdx = new int[num_springs * 4]; {
		for (int i = 0; i < num_springs * 4; i++) {
			sprIdx[i] = spring_connectivity[i];
		}
	}


	//tri neighbours
	///BW addition -> first get triangle neighbours for each node
	tri_neighbours = new int[num_nodes *maxFreedom];
	tri_neighbours_minor = new int[num_nodes *maxFreedom];
	minor_node_neighbour = new int[num_nodes *maxFreedom];
	get_tet_hash();
	find_tet_neighbours();

	return;

}



//ported from Ming Zhu / Osay code to allow for GPU compatible memory types
//not sure what motivation behind some tasks is

void lagrangian_object::Mesh_Process() {

	std::cout << "\tProcessing mesh file..." << std::endl;



	//step 1: find nodFreedom according to the number of triangles 
	//which share the node.

	//understood
	int tempNodFreedom[num_nodes];
	for (int i = 0; i < num_nodes; i++)
		tempNodFreedom[i] = 0;
	for (int i = 0; i < num_tets; i++) {
		for (int j = 0; j < 3; j++) {
			tempNodFreedom[tet_connectivity[i * 3 + j]]++;
		}
	}

	//step 2: find the max and min nodfreedom of the spring network. //understood
	maxFreedom = tempNodFreedom[0];
	minFreedom = tempNodFreedom[0];
	for (int i = 1; i < num_nodes; i++) {
		maxFreedom = (maxFreedom > tempNodFreedom[i]) ? maxFreedom : tempNodFreedom[i];
		minFreedom = (minFreedom < tempNodFreedom[i]) ? minFreedom : tempNodFreedom[i];
	}


	//step 3: find index of triangles around each node. //understood
	int temp_nodTriIdx[num_nodes][maxFreedom];
	int nodTriCount[num_nodes];
	for (int i = 0; i < num_nodes; i++)
		nodTriCount[i] = 0;
	for (int i = 0; i < num_tets; i++) {
		for (int j = 0; j < 3; j++) {
			const int k = tet_connectivity[i * 3 + j];

			temp_nodTriIdx[k][nodTriCount[k]++] = i;
		}
	}

	//step 4: find index of nodes around each node.
	int nodTriIdx[num_nodes][maxFreedom];

	int** tempNodIdx = new int*[num_nodes];

	for (int i = 0; i < num_nodes; i++) {

		const int i_dfm = tempNodFreedom[i];

		//std::cout << i+1 << "\t" << i_dfm << "\t";

		//extract all the nodIdx around the nod_i;
		int nodIdx_local[i_dfm * 2];
		for (int j = 0; j < i_dfm; j++) {

			const int triIdx_local = temp_nodTriIdx[i][j];

			//std::cout << triIdx_local << "\t";

			for (int k = 0; k < 3; k++) {
				if (tet_connectivity[triIdx_local * 3 + k] == i) {

					nodIdx_local[j * 2] = tet_connectivity[triIdx_local * 3 + boundaryIdx(k + 1, 3)];
					nodIdx_local[j * 2 + 1] = tet_connectivity[triIdx_local * 3 + boundaryIdx(k + 2, 3)];

					break;
				}

			}

		}

		/// nodIdx_local holds all neighbours from each of the neighbouring triangles

		tempNodIdx[i] = new int[i_dfm + 1];

		int tempNodIdx_count = -1;
		tempNodIdx[i][++tempNodIdx_count] = nodIdx_local[0];
		tempNodIdx[i][++tempNodIdx_count] = nodIdx_local[1];

		bool triIdx_local_register[i_dfm];
		triIdx_local_register[0] = false;
		nodTriIdx[i][0] = temp_nodTriIdx[i][0];


		/// guess that the rest of this is to do "edge" stuff
		/// not relevent to RBC calcs

		for (int j = 1; j < i_dfm; j++)
			triIdx_local_register[j] = true;

		for (int j = 1; j < i_dfm; j++) {
			bool edge = true;

			for (int k = 1; k < i_dfm; k++) {
				if (triIdx_local_register[k]) {
					if (nodIdx_local[k * 2] == tempNodIdx[i][tempNodIdx_count]) {

						nodTriIdx[i][tempNodIdx_count] = temp_nodTriIdx[i][k];
						tempNodIdx[i][++tempNodIdx_count] = nodIdx_local[k * 2 + 1];

						triIdx_local_register[k] = false;
						edge = false;
						break;
					}
					else if (nodIdx_local[k * 2 + 1] == tempNodIdx[i][tempNodIdx_count]) {

						nodTriIdx[i][tempNodIdx_count] = temp_nodTriIdx[i][k];
						tempNodIdx[i][++tempNodIdx_count] = nodIdx_local[k * 2];
						triIdx_local_register[k] = false;
						edge = false;
						break;
					}
				}
			}

			if (edge) {

				//std::cout << i+1 << "\t" << i_dfm << "\t";




				const int nod_remain = i_dfm - j;
				for (int k = i_dfm; k > nod_remain; k--) {
					nodTriIdx[i][k - 1] = nodTriIdx[i][k - 1 - nod_remain];
					tempNodIdx[i][k] = tempNodIdx[i][k - nod_remain];
				}
				tempNodIdx[i][nod_remain] = tempNodIdx[i][0];

				tempNodIdx_count = nod_remain;

				for (int k = 0; k < nod_remain; k++) {

					for (int w = 1; w < i_dfm; w++) {
						if (triIdx_local_register[w]) {
							if (nodIdx_local[w * 2] == tempNodIdx[i][tempNodIdx_count]) {

								tempNodIdx[i][--tempNodIdx_count] = nodIdx_local[w * 2 + 1];
								nodTriIdx[i][tempNodIdx_count] = temp_nodTriIdx[i][w];
								triIdx_local_register[w] = false;
								break;
							}
							else if (nodIdx_local[w * 2 + 1] == tempNodIdx[i][tempNodIdx_count]) {

								tempNodIdx[i][--tempNodIdx_count] = nodIdx_local[w * 2];
								nodTriIdx[i][tempNodIdx_count] = temp_nodTriIdx[i][w];
								triIdx_local_register[w] = false;
								break;
							}
						}
					}
				}



				j = i_dfm; //secure.
				break;
			}
		}


	}

	//step 5: sort out the node index and triangle index.
	//and distinguish the nodes at edge. 
	//importantly, update the nodFreedom according to the number of springs 
	//which share a node.

	/// BW: sorting goes on here

	bool temp_nodEdge[num_nodes];

	int nodIdx_sequence[num_nodes][2];
	for (int i = 0; i < num_nodes; i++)
		nodIdx_sequence[i][0] = nodIdx_sequence[i][1] = -1;
	int nodIdx_sequence_count = 0;

	int triIdx_sequence[num_tets][2];
	for (int i = 0; i < num_tets; i++)
		triIdx_sequence[i][0] = triIdx_sequence[i][1] = -1;
	int triIdx_sequence_count = 1;
	triIdx_sequence[0][0] = nodTriIdx[0][0];
	triIdx_sequence[nodTriIdx[0][0]][1] = 0;

	for (int i = 0; i < num_tets; i++) {

		for (int j = 0; j < 3; j++) {

			const int i_triIdx = triIdx_sequence[i][0];

			const int j_nodIdx = tet_connectivity[i_triIdx * 3 + j];

			const int j_dfm = tempNodFreedom[j_nodIdx];

			//check if the node is registered.
			if (nodIdx_sequence[j_nodIdx][1] == -1) {
				//register.
				nodIdx_sequence[nodIdx_sequence_count][0] = j_nodIdx;
				nodIdx_sequence[j_nodIdx][1] = nodIdx_sequence_count;

				nodIdx_sequence_count++;
				//check the order of nodes around the node j_nodIdx.
				bool ordered = false;
				for (int k = 0; k < j_dfm; k++) {

					if ((tempNodIdx[j_nodIdx][k] == tet_connectivity[i_triIdx * 3 + boundaryIdx(j + 1, 3)])   //boundary idx is mod 3 function
						&& (tempNodIdx[j_nodIdx][k + 1] == tet_connectivity[i_triIdx * 3 + boundaryIdx(j + 2, 3)])) {


						ordered = true;
						break;
					}
				}

				if (ordered == false) {

					std::cout << "reversing..." << std::endl;;
					

					for (int w = 0; w < (j_dfm + 1) / 2; w++) {
						int tw = nodTriIdx[j_nodIdx][w];
						tempNodIdx[j_nodIdx][w] = tempNodIdx[j_nodIdx][j_dfm - 1 - w];
						tempNodIdx[j_nodIdx][j_dfm - 1 - w] = tw;

						tw = tempNodIdx[j_nodIdx][w];
						tempNodIdx[j_nodIdx][w] = tempNodIdx[j_nodIdx][j_dfm - w];
						tempNodIdx[j_nodIdx][j_dfm - w] = tw;
					}
				}

				//check if triangles around the node j_nodIndex are registered.
				for (int k = 0; k < j_dfm; k++) {

					const int k_triIdx = nodTriIdx[j_nodIdx][k];

					if (triIdx_sequence[k_triIdx][1] == -1) {


						triIdx_sequence[triIdx_sequence_count][0] = k_triIdx;
						triIdx_sequence[k_triIdx][1] = triIdx_sequence_count;

						triIdx_sequence_count++;

					}
				}


				if (tempNodIdx[j_nodIdx][j_dfm] != tempNodIdx[j_nodIdx][0]) {
					num_springs += (++tempNodFreedom[j_nodIdx]);
					temp_nodEdge[j_nodIdx] = true;

					//std::cout << "edge" << j_nodIdx
				}
				else {
					num_springs += (tempNodFreedom[j_nodIdx]);
					temp_nodEdge[j_nodIdx] = false;
				}

			}
		}
	}


	//step 6: AGAIN update the max and min nodfreedom of the spring network.
	maxFreedom = tempNodFreedom[0];
	minFreedom = tempNodFreedom[0];
	for (int i = 1; i < num_nodes; i++) {
		maxFreedom = (maxFreedom > tempNodFreedom[i]) ? maxFreedom : tempNodFreedom[i];
		minFreedom = (minFreedom < tempNodFreedom[i]) ? minFreedom : tempNodFreedom[i];
	}

	//step 7: finalise the nodIdx, nodFreedom, nodCoordinate, triIdx.
	nodIdx = new int[num_nodes* maxFreedom];
	if (nodIdx == NULL) exit(1);
	nodEdge = new bool[num_nodes];
	if (nodEdge == NULL) exit(1);
	sorted_nodFreedom = new int[num_nodes];
	if (sorted_nodFreedom == NULL) exit(1);
	sorted_coordinate = new double[num_nodes * 3];

	for (int i = 0; i < num_nodes; i++) {
		const int i_nodIdx = nodIdx_sequence[i][0];

		const int i_dfm = tempNodFreedom[i_nodIdx];

		for (int j = 0; j < i_dfm; j++) {

			nodIdx[i*maxFreedom + j] = nodIdx_sequence[tempNodIdx[i_nodIdx][j]][1];

		}
		sorted_nodFreedom[i] = i_dfm;

		sorted_coordinate[i * 3] = node_x[i_nodIdx];
		sorted_coordinate[i * 3 + 1] = node_y[i_nodIdx];
		sorted_coordinate[i * 3 + 2] = node_z[i_nodIdx];

		nodEdge[i] = temp_nodEdge[i_nodIdx];
	}


	sorted_triIdx = new int[num_tets * 3];
	for (int i = 0; i < num_tets; i++) {
		const int i_triIdx = triIdx_sequence[i][0];

		for (int j = 0; j < 3; j++) {

			sorted_triIdx[i * 3 + j] = nodIdx_sequence[tet_connectivity[i_triIdx * 3 + j]][1];

		}
	}



	//step 9: finalise the sprIdx.
	if (num_springs % 2 == 0) {
		num_springs /= 2;
	}
	else {
		std::cout << "Error code: 1050." << std::endl;
		getchar();
	}

	sprIdx = new int[num_springs * 4];
	sprEdge = new bool[num_springs];

	for (int i = 0; i < num_springs; i++)
		sprEdge[i] = false;

	int sprN_counter = 0;
	for (int i = 0; i < num_nodes; i++) {

		const int i_dfm = sorted_nodFreedom[i];


		for (int j = 0; j < i_dfm; j++) {

			if (i < nodIdx[i * maxFreedom + j]) {
				sprIdx[sprN_counter *4 + 0] = i;
				sprIdx[sprN_counter * 4 + 1] = nodIdx[i*maxFreedom + boundaryIdx(j - 1, i_dfm)];
				sprIdx[sprN_counter * 4 + 2] = nodIdx[i * maxFreedom + j];
				sprIdx[sprN_counter * 4 + 3] = nodIdx[i * maxFreedom + boundaryIdx(j + 1, i_dfm)];

				if ((nodEdge[i] == true) && (j == 0)) {
					sprIdx[sprN_counter * 4 + 1] = sprIdx[sprN_counter * 4 + 3];
					sprEdge[sprN_counter * 4 + 0] = true;
				}
				else if ((nodEdge[i] == true) && (j + 1 == i_dfm)) {
					sprIdx[sprN_counter * 4 + 3] = sprIdx[sprN_counter * 4 + 1];
					sprEdge[sprN_counter * 4 + 0] = true;
				}

				sprN_counter++;
			}
		}
	}

	//step 10: determine the mesh is closed.
	closed_mesh = true;
	for (int i = 0; i < num_nodes; i++) {
		if (nodEdge[i] == true) {
			closed_mesh = false;
			break;
		}
	}


	///BW addition -> first get triangle neighbours for each node
	tri_neighbours = new int[num_nodes *maxFreedom];
	tri_neighbours_minor = new int[num_nodes *maxFreedom];
	minor_node_neighbour = new int[num_nodes *maxFreedom];
	get_tet_hash();
	find_tet_neighbours();

	



	for (int i = 0; i < num_nodes; i++)
		delete[] tempNodIdx[i];
	delete[] tempNodIdx;

	return;
}


void lagrangian_object::find_tet_neighbours() {

	int node_dof;

	int neighbour;

	std::vector <int> temp;
	long long int hash_int;
	int t, of_face;
	int factor;

	factor = ceil(log10(num_nodes));

	int w, x, y, z;

	//create hash from Node neighbours array

	for (int i = 0; i < num_nodes; i++) {
		node_dof = sorted_nodFreedom[i];

		//do for first ndof -1 triangle

		for (int j = 0; j < node_dof ; j++) {
			temp.clear();

			temp.push_back(i);
			temp.push_back(nodIdx[i*maxFreedom + j] );
			temp.push_back(nodIdx[i*maxFreedom + ((j + 1) % node_dof + node_dof) % node_dof]);

			std::sort(temp.begin(), temp.end());
			hash_int = 0;
			t = 0;

			for (std::vector<int>::iterator it = temp.begin(); it != temp.end(); ++it) {

				hash_int = hash_int + *it * pow(10, factor * t);
				t++;
			}


			std::vector<long long int>::iterator it;

			int p;
			
			//find hash_int in tet_hash
			
				it = std::find(tet_hash.begin(), tet_hash.end(), hash_int);
				if (it != tet_hash.end()) {
					p = std::distance(tet_hash.begin(), it);
					tri_neighbours[i*maxFreedom + j] = p;
				}


			//now to find minor index

			//first get third node of minor 

			x = nodIdx[i*maxFreedom + j];
			y = nodIdx[i*maxFreedom + ((j + 1) % node_dof + node_dof) % node_dof];

			for (int k = 0; k < sorted_nodFreedom[x]; k++) {
				if (nodIdx[x*maxFreedom + k] == y) {
					z = k;
					break;  // z is y node
				}
			}
			
			if (nodIdx[x*maxFreedom + ((z - 1) % sorted_nodFreedom[x] + sorted_nodFreedom[x]) % sorted_nodFreedom[x]] == i) {
				w = nodIdx[x*maxFreedom + ((z + 1) % sorted_nodFreedom[x] + sorted_nodFreedom[x]) % sorted_nodFreedom[x]];
			}
			else {
				w = nodIdx[x*maxFreedom + ((z - 1) % sorted_nodFreedom[x] + sorted_nodFreedom[x]) % sorted_nodFreedom[x]];
			}

			minor_node_neighbour[i*maxFreedom + j] = w;

			temp.clear();

			temp.push_back(x);
			temp.push_back(y);
			temp.push_back(w);

			std::sort(temp.begin(), temp.end());
			hash_int = 0;
			t = 0;

			for (std::vector<int>::iterator it = temp.begin(); it != temp.end(); ++it) {

				hash_int = hash_int + *it * pow(10, factor * t);
				t++;
			}


			//find hash_int in tet_hash

			it = std::find(tet_hash.begin(), tet_hash.end(), hash_int);
			if (it != tet_hash.end()) {
				p = std::distance(tet_hash.begin(), it);
				tri_neighbours_minor[i*maxFreedom + j] = p;
			}
			
		}

	}

}



void lagrangian_object::get_tet_hash() {

	int neighbour;

	std::vector <int> temp;
	long long int hash_int;
	int t, of_face;
	int factor;

	factor = ceil(log10(num_nodes));

	for (int i = 0; i < num_tets; i++) {

		temp.clear();

		temp.push_back(sorted_triIdx[i*3] );
		temp.push_back(sorted_triIdx[i*3 +1]);
		temp.push_back(sorted_triIdx[i *3 +2]);

		std::sort(temp.begin(), temp.end());
		hash_int = 0;
		t = 0;

		for (std::vector<int>::iterator it = temp.begin(); it != temp.end(); ++it) {

			hash_int = hash_int + *it * pow(10, factor * t);
			t++;
		}

		tet_hash.push_back(hash_int);
		// then identify the ghost cell parameters needed

	}

}
void lagrangian_object::solution_read(global_variables &globals) {

	std::cout << "\tImporting Existing Solution..." << std::endl;

	//just neeed x,y,z coordinates . Force and Area calculated in solution

	std::ifstream inputFile;
	std::stringstream tt;
	std::string inputFileName;
	tt << ceil(globals.simulation_length);
	inputFileName = globals.output_file + "/" + name + ".dat";

	inputFile.open(inputFileName.c_str());
	std::string line;
	std::string dummy;

	if (inputFile.is_open()) {
		getline(inputFile, line);
		getline(inputFile, line);
		getline(inputFile, line);
		getline(inputFile, line);

		for (int i = 0; i < num_nodes; i++) {

			getline(inputFile, line);

			std::stringstream   linestream(line);
			/*	std::string         data;
				int                 val1;
				int                 val2;*/

			linestream >> node_x[i];
			linestream >> node_y[i];
			linestream >> node_z[i];
			linestream >> node_force_x[i];
			linestream >> node_force_y[i];
			linestream >> node_force_z[i];
			linestream >> node_vel_x[i];
			linestream >> node_vel_y[i];
			linestream >> node_vel_z[i];
			linestream >> dummy;

		}


		std::cout << "\t .... Solution Imported" << std::endl;



	}
	else {
		std::cout << "Failed to open the mesh solution file. Press ENTER to exit!" << std::endl;
		getchar();
		exit(EXIT_FAILURE);
	}
	inputFile.close();
}


void lagrangian_object::read_initial_mesh(std::string inputFileName) {

	std::ifstream inputFile;
	inputFile.open(inputFileName.c_str());
	if (inputFile.is_open()) {

		//std::cout << "Successfully opened the NODES input file! " << std::endl;

		//for nodes.
		inputFile >> num_nodes;
		inputFile >> num_tets;

		//std::cout << num_nodes << " " << triN << std::endl;

		//update x_y_z nodes
		node_x = new double[num_nodes];
		if (node_x == NULL) exit(1);
		node_y = new double[num_nodes];
		if (node_y == NULL) exit(1);
		node_z = new double[num_nodes];
		if (node_z == NULL) exit(1);
		node_x_ref = new double[num_nodes];
		if (node_x_ref == NULL) exit(1);
		node_y_ref = new double[num_nodes];
		if (node_y_ref == NULL) exit(1);
		node_z_ref = new double[num_nodes];
		if (node_z_ref == NULL) exit(1);

		//import xyz nodes
		for (int i = 0; i < num_nodes; i++) {

			inputFile >> node_x[i];
			inputFile >> node_y[i];
			inputFile >> node_z[i];
			node_x_ref[i] = node_x[i] + centre_x;
			node_y_ref[i] = node_y[i] + centre_y;
			node_z_ref[i] = node_z[i] + centre_z;
			node_x[i] = node_x_ref[i];
			node_y[i] = node_y_ref[i];
			node_z[i] = node_z_ref[i];

		}
		
	}
	else {
		std::cout << "Failed to open the mesh solution file. Press ENTER to exit!" << std::endl;
		getchar();
		exit(EXIT_FAILURE);
	}
	inputFile.close();
}


void lagrangian_object::read_initial_mesh_MZ(std::string inputFileName) {

	std::ifstream inputFile;
	inputFile.open(inputFileName.c_str());
	if (inputFile.is_open()) {

		//std::cout << "Successfully opened the NODES input file! " << std::endl;
//for nodes.
		inputFile >> num_nodes;
		inputFile >> num_springs;
		inputFile >> num_tets;

		//std::cout << num_nodes << " " << triN << std::endl;

		//update x_y_z nodes
		
		//import xyz nodes
		for (int i = 0; i < num_nodes; i++) {

			inputFile >> node_x[i];
			inputFile >> node_y[i];
			inputFile >> node_z[i];
			node_x_ref[i] = node_x[i] + centre_x;
			node_y_ref[i] = node_y[i] + centre_y;
			node_z_ref[i] = node_z[i] + centre_z;
			node_x[i] = node_x_ref[i];
			node_y[i] = node_y_ref[i];
			node_z[i] = node_z_ref[i];


		}

	}
	else {
		std::cout << "Failed to open the mesh solution file. Press ENTER to exit!" << std::endl;
		getchar();
		exit(EXIT_FAILURE);
	}
	inputFile.close();
}

