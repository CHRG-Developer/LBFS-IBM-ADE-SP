#pragma once
#include <string>
#include <vector>
#include <algorithm>
#include "global_variables.h"

using namespace std;
class lagrangian_object
{
public:
	lagrangian_object();
	~lagrangian_object();

	std::string name, mesh_file, stress_free_mesh, curvature_mesh;
	int num_nodes, num_springs, num_tets;
	int min_index, max_index;
	int type; // 1 = rigid_cylinder
	double radius, depth;
	double centre_x, centre_y, centre_z;
	double stiffness;
	double total_volume, total_volume_0;
	double total_area, total_area_0;
	double total_MAD, total_MAD_0;
	double boundary_force;
	double point_mass;
	
	double * node_x, *node_y, *node_z;
	double * node_x_ref, *node_y_ref, *node_z_ref;
	double * node_vel_x, *node_vel_y, *node_vel_z;
	double * node_force_x, *node_force_y, *node_force_z;
	int depth_nodes, radial_nodes;

	int * tet_connectivity, *spring_connectivity;
	int * num_node_neighbours , * minor_node_neighbour;
	int * nodIdx, *sorted_triIdx, *tri_neighbours, *tri_neighbours_minor;  // index of neighbouring nodes and triangles for a node
	int * spring_neighbours;
	int * sprIdx;
	bool *sprEdge, *nodEdge; // is spring or node an edge  - don't know what an edge is
	double * external_loading_factor;
	bool periodic = false;
	bool sphere_reference = true;
	bool mz_importer = false;
	double spontaneous_angle;
	double volume_ratio = 1.0;
	double reference_length;
	bool pivkin_length = false;

	int maxFreedom, minFreedom;
	std::vector <std::vector <int> > node_neighbours;
	vector<double>	nodCoordinate;
	vector<double>	nodCoordinate_initial;
	std::vector <long long int>  tet_hash;

	int * sorted_nodFreedom, *nodFreedom;
	double * sorted_coordinate;

	void initialise(double PI);

	void populate_node_reference_displacement(double PI);
	void import_network_mesh();
	void mesh_read(std::string inputFileName);
	void osay_mesh_read(std::string inputFileName);
	void read_initial_mesh(std::string inputFileName);
	void solution_read(global_variables &globals);
	void pre_process_spring_network();
	void Mesh_Process();
	void add_volume(double* tet_volume_block, bool first_timestep, int n_blocks);
	void add_area(double* vol, bool first_timestep, int n_blocks);
	void add_MAD(double* MAD, bool first_timestep, int n_blocks);
	void get_tet_hash();
	void find_tet_neighbours();
	void apply_optical_tweezers_force();
	void compatibility();
	void read_initial_mesh_MZ(std::string inputFileName);

	double volume_modulus;
	double local_area_modulus;
	double global_area_modulus;
	double shear_modulus;
	double global_bending_modulus;
	double local_bending_modulus;
	double membrane_thickness =1.0;
	double viscosity;
	double wlc_ratio;
	double internal_viscosity_ratio =5.0;
	double curvature_ratio;


	// Ming Zhu function - recursive to find 

	// BW - finds x mod Idx  that's between 0 and Idx 
	int boundaryIdx(int x, int idx) {
		if (x >= idx)
			return boundaryIdx(x - idx, idx);
		else if (x < 0)
			return boundaryIdx(x + idx, idx);
		else
			return x;
	}


	
	


private:
	bool			closed_mesh;
};

struct target_less
{
	template<class It>
	bool operator()(It const &a, It const &b) const { return *a < *b; }
};
struct target_equal
{
	template<class It>
	bool operator()(It const &a, It const &b) const { return *a == *b; }
};
template<class It> It uniquify(It begin, It const end)
{
	std::vector<It> v;
	v.reserve(static_cast<size_t>(std::distance(begin, end)));
	for (It i = begin; i != end; ++i)
	{
		v.push_back(i);
	}
	std::sort(v.begin(), v.end(), target_less());
	v.erase(std::unique(v.begin(), v.end(), target_equal()), v.end());
	std::sort(v.begin(), v.end());
	size_t j = 0;
	for (It i = begin; i != end && j != v.size(); ++i)
	{
		if (i == v[j])
		{
			using std::iter_swap; iter_swap(i, begin);
			++j;
			++begin;
		}
	}
	return begin;
}