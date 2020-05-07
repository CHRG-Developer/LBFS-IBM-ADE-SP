#include "preprocessor.h"
#include "global_variables.h"
#include <fstream>
#include "tinyxml2.h"
#include <sstream>
#include "initial_conditions.h"
#include "quad_bcs_plus.h"
#include <vector>
#include <iostream>
#include <cstdio>
#include <algorithm>
#include "lagrangian_object.h"

using namespace tinyxml2;

preprocessor::preprocessor()
{
    //ctor

}

preprocessor::~preprocessor()
{
    //dtor
}

void preprocessor::initialise_program_variables(char* xmltest, global_variables &globals, domain_geometry &geometry,
                                                initial_conditions &initial_conds, quad_bcs_plus &bcs,
                                                unstructured_bcs &u_bcs, vector<lagrangian_object> &object_vec) {

    XMLDocument xmlDoc;

    XMLError eResult = xmlDoc.LoadFile(xmltest);


    parse_global_variables(xmlDoc,  globals);
    parse_geometry_variables(xmlDoc,  geometry);
    parse_initial_conditions(xmlDoc, initial_conds);
    parse_boundary_conditions(xmlDoc, bcs);
    //unstructured mesh types
    if (globals.mesh_type > 2){
        parse_uns_boundary_conditions(xmlDoc, u_bcs,globals);
    }
	parse_lagrangian_objects(xmlDoc, object_vec, globals);

    mach_number_factor(globals,bcs,initial_conds,geometry);
    geometry.scale_geometries(globals.scale);
    globals.initialise(geometry,initial_conds);



}
void preprocessor::mach_number_factor( global_variables &globals,quad_bcs_plus &bcs,
        initial_conditions &initials,domain_geometry &geometry ){


    double factor = globals.max_velocity/globals.scale ;
    //double factor = globals.max_velocity;

    bcs.e_u = bcs.e_u * factor;
    bcs.w_u = bcs.w_u * factor;
    bcs.n_u = bcs.n_u * factor;
    bcs.s_u = bcs.s_u *factor;

    bcs.e_v = bcs.e_v * factor;
    bcs.w_v = bcs.w_v * factor;
    bcs.n_v = bcs.n_v * factor;
    bcs.s_v = bcs.s_v *factor;

    initials.vel_gradient.x = initials.vel_gradient.x * factor ;
    initials.vel_gradient.y = initials.vel_gradient.y * factor;
    initials.vel_gradient.z = initials.vel_gradient.z * factor;
    initials.vel_origin_mag.x = initials.vel_origin_mag.x * factor;
    initials.vel_origin_mag.y = initials.vel_origin_mag.y * factor;
    initials.vel_origin_mag.z = initials.vel_origin_mag.z * factor;

    initials.velocity.x = initials.velocity.x * factor;
    initials.velocity.y = initials.velocity.y * factor;
    initials.velocity.z = initials.velocity.z * factor;



}



void preprocessor::parse_lagrangian_objects(XMLDocument &xmlDoc, std::vector<lagrangian_object> &object_vec, global_variables &globals) {

	XMLError eResult;

	XMLNode * pRoot = xmlDoc.FirstChild();
	const char * output;
	std::string temp;
	std::ostringstream os;
	int temp_type;
	if (pRoot == nullptr) return;


	double iOutListValue;

	XMLElement * pElement = pRoot->FirstChildElement("lagrangian_objects");
	if (pElement == nullptr) return;

	XMLElement * pListElement = pElement->FirstChildElement("object");
	XMLElement * pbcElement;
	lagrangian_object temp_obj;

	while (pListElement != nullptr)
	{

		pbcElement = pListElement->FirstChildElement("name");
		if (pbcElement == nullptr) return;
		os << pbcElement->GetText();
		temp_obj.name = os.str();

		os.str(std::string());
		os.clear();


		pbcElement = pListElement->FirstChildElement("num_nodes");
		eResult = pbcElement->QueryDoubleText(&iOutListValue);
		temp_obj.num_nodes = iOutListValue;

		pbcElement = pListElement->FirstChildElement("radius");
		eResult = pbcElement->QueryDoubleText(&iOutListValue);
		temp_obj.radius = iOutListValue;

		pbcElement = pListElement->FirstChildElement("depth");
		eResult = pbcElement->QueryDoubleText(&iOutListValue);
		temp_obj.depth = iOutListValue;

		pbcElement = pListElement->FirstChildElement("stiffness");
		eResult = pbcElement->QueryDoubleText(&iOutListValue);
		temp_obj.stiffness = iOutListValue;

		pbcElement = pListElement->FirstChildElement("centre_x");
		eResult = pbcElement->QueryDoubleText(&iOutListValue);
		temp_obj.centre_x = iOutListValue;

		pbcElement = pListElement->FirstChildElement("centre_y");
		eResult = pbcElement->QueryDoubleText(&iOutListValue);
		temp_obj.centre_y = iOutListValue;

		pbcElement = pListElement->FirstChildElement("centre_z");
		eResult = pbcElement->QueryDoubleText(&iOutListValue);
		temp_obj.centre_z = iOutListValue;

		pbcElement = pListElement->FirstChildElement("type");
		os << pbcElement->GetText();

		temp = os.str();
		if (temp.compare("rigid_cylinder") == 0) {
			temp_obj.type = 1;
		}else if (temp.compare("spring_network") == 0) {
			temp_obj.type = 2;
		}
		else  {
			temp_obj.type = 0;
		}

		os.clear();
		os.str(std::string());

		//get mesh file of spring network

		pbcElement = pListElement->FirstChildElement("mesh_file");
		if (pbcElement == nullptr) {

		}else {
			os << pbcElement->GetText();
			temp_obj.mesh_file = os.str();


		}

		pbcElement = pListElement->FirstChildElement("local_area_modulus");
		if (pbcElement == nullptr) {

		}else {
			eResult = pbcElement->QueryDoubleText(&iOutListValue);
			temp_obj.local_area_modulus = iOutListValue;
		}

		pbcElement = pListElement->FirstChildElement("global_area_modulus");
		if (pbcElement == nullptr) {

		}
		else {
			eResult = pbcElement->QueryDoubleText(&iOutListValue);
			temp_obj.global_area_modulus = iOutListValue;
		}

		pbcElement = pListElement->FirstChildElement("volume_modulus");
		if (pbcElement == nullptr) {

		}
		else {
			eResult = pbcElement->QueryDoubleText(&iOutListValue);
			temp_obj.volume_modulus = iOutListValue;
		}

		pbcElement = pListElement->FirstChildElement("shear_modulus");
		if (pbcElement == nullptr) {

		}
		else {
			eResult = pbcElement->QueryDoubleText(&iOutListValue);
			temp_obj.shear_modulus = iOutListValue;
		}

		pbcElement = pListElement->FirstChildElement("local_bending_modulus");
		if (pbcElement == nullptr) {

		}
		else {
			eResult = pbcElement->QueryDoubleText(&iOutListValue);
			temp_obj.local_bending_modulus = iOutListValue;
		}

		pbcElement = pListElement->FirstChildElement("global_bending_modulus");
		if (pbcElement == nullptr) {

		}
		else {
			eResult = pbcElement->QueryDoubleText(&iOutListValue);
			temp_obj.global_bending_modulus = iOutListValue;
		}

		pbcElement = pListElement->FirstChildElement("membrane_thickness");
		if (pbcElement == nullptr) {

		}
		else {
			eResult = pbcElement->QueryDoubleText(&iOutListValue);
			temp_obj.membrane_thickness = iOutListValue;
		}

		pbcElement = pListElement->FirstChildElement("viscosity");
		if (pbcElement == nullptr) {

		}
		else {
			eResult = pbcElement->QueryDoubleText(&iOutListValue);
			temp_obj.viscosity = iOutListValue;
		}

		pbcElement = pListElement->FirstChildElement("wlc_ratio");
		if (pbcElement == nullptr) {

		}
		else {
			eResult = pbcElement->QueryDoubleText(&iOutListValue);
			temp_obj.wlc_ratio = iOutListValue;
		}

		pbcElement = pListElement->FirstChildElement("boundary_force");
		if (pbcElement == nullptr) {

		}
		else {
			eResult = pbcElement->QueryDoubleText(&iOutListValue);
			temp_obj.boundary_force = iOutListValue;
		}

		pbcElement = pListElement->FirstChildElement("point_mass");
		if (pbcElement == nullptr) {

		}
		else {
			eResult = pbcElement->QueryDoubleText(&iOutListValue);
			temp_obj.point_mass = iOutListValue;
		}
		pbcElement = pListElement->FirstChildElement("internal_viscosity_ratio");
		if (pbcElement == nullptr) {

		}
		else {
			eResult = pbcElement->QueryDoubleText(&iOutListValue);
			temp_obj.internal_viscosity_ratio = iOutListValue;
		}

		os.clear();
		os.str(std::string());

		pbcElement = pListElement->FirstChildElement("periodic");
		if (pbcElement == nullptr) {

		}
		else {
			os << pbcElement->GetText();
		}
	
		temp = os.str();
		if (temp.compare("true") == 0) {
			temp_obj.periodic = true;
		}
		else {
			temp_obj.periodic = false;
		}

		pbcElement = pListElement->FirstChildElement("volume_ratio");
		if (pbcElement == nullptr) {

		}
		else {
			eResult = pbcElement->QueryDoubleText(&iOutListValue);
			temp_obj.volume_ratio = iOutListValue;
		}

		os.clear();
		os.str(std::string());

		//get mesh file of spring network

		pbcElement = pListElement->FirstChildElement("stress_free_mesh");
		if (pbcElement == nullptr) {

		}
		else {
			os << pbcElement->GetText();
			temp_obj.stress_free_mesh = os.str();
		}

		os.clear();
		os.str(std::string());

		//get mesh file of spring network

		pbcElement = pListElement->FirstChildElement("curvature_mesh");
		if (pbcElement == nullptr) {

		}
		else {
			os << pbcElement->GetText();
			temp_obj.curvature_mesh = os.str();
		}

		os.clear();
		os.str(std::string());

		pbcElement = pListElement->FirstChildElement("sphere_reference");
		if (pbcElement == nullptr) {

		}
		else {
			os << pbcElement->GetText();
		}

		temp = os.str();
		if (temp.compare("true") == 0) {
			temp_obj.periodic = true;
		}
		else {
			temp_obj.periodic = false;
		}


		pbcElement = pListElement->FirstChildElement("spring_length");
		if (pbcElement == nullptr) {

		}
		else {
			eResult = pbcElement->QueryDoubleText(&iOutListValue);
			if (iOutListValue < 0) {
				temp_obj.reference_length = iOutListValue;
				temp_obj.pivkin_length = false;
			}
			else {
				temp_obj.reference_length = iOutListValue;
				temp_obj.pivkin_length = true;
			}
			
		}

		os.clear();
		os.str(std::string());

		pbcElement = pListElement->FirstChildElement("mz_importer");
		if (pbcElement == nullptr) {

		}
		else {
			os << pbcElement->GetText();
		}

		temp = os.str();
		if (temp.compare("true") == 0) {
			temp_obj.mz_importer = true;
		}
		else {
			temp_obj.mz_importer = false;
		}


		pbcElement = pListElement->FirstChildElement("curvature_ratio");
		if (pbcElement == nullptr) {

		}
		else {
			eResult = pbcElement->QueryDoubleText(&iOutListValue);
			temp_obj.curvature_ratio = iOutListValue;
		}


		//temp_obj.initialise(globals.PI);

		object_vec.push_back(temp_obj);

		//loop to next lagnrangian object condition
		pListElement = pListElement->NextSiblingElement("object");

	}

}


void preprocessor::parse_uns_boundary_conditions(XMLDocument &xmlDoc, unstructured_bcs &bcs,
            global_variables &globals){

    double factor = globals.max_velocity/globals.scale ;

    XMLError eResult;

    XMLNode * pRoot = xmlDoc.FirstChild();
    const char * output;
    std::string temp;
   std::ostringstream os;
   int temp_type;
    if (pRoot == nullptr) return;


        double iOutListValue;

        XMLElement * pElement = pRoot->FirstChildElement("uns_boundary_conditions");
        if (pElement == nullptr) return ;

        XMLElement * pListElement = pElement->FirstChildElement("bc");
         XMLElement * pbcElement;

         while (pListElement != nullptr)
            {
                pbcElement = pListElement->FirstChildElement("name");
                if (pbcElement == nullptr) return ;
                os << pbcElement->GetText();
                bcs.name.push_back(os.str());
                std::cout << "name of BC " << os.str() << std::endl;
                os.str(std::string());
                os.clear();

                pbcElement = pListElement->FirstChildElement("rho");
                eResult = pbcElement->QueryDoubleText(&iOutListValue);
                bcs.rho.push_back(iOutListValue );

                pbcElement = pListElement->FirstChildElement("u");
                eResult = pbcElement->QueryDoubleText(&iOutListValue);
                bcs.u.push_back(iOutListValue*factor);

                pbcElement = pListElement->FirstChildElement("v");
                eResult = pbcElement->QueryDoubleText(&iOutListValue);
                bcs.v.push_back(iOutListValue*factor);

                pbcElement = pListElement->FirstChildElement("w");
                eResult = pbcElement->QueryDoubleText(&iOutListValue);
                bcs.w.push_back(iOutListValue*factor);


                pbcElement = pListElement->FirstChildElement("vel_type");
                os << pbcElement->GetText();

                temp = os.str();
                if (temp.compare("dirichlet") == 0){
                    temp_type= 1;
                }else if ( temp.compare("neumann") == 0){
                    temp_type = 2;
                }else if ( temp.compare("periodic") == 0){
                    temp_type = 3;
                }else if ( temp.compare("parabolic-N") == 0){
                    temp_type = 4;
                }else if ( temp.compare("parabolic-W") == 0){
                    temp_type= 5;
                }else if ( temp.compare("inlet") == 0){
                    temp_type = 6;
                }else if ( temp.compare("wall") == 0){
                    temp_type = 7;
                }else if ( temp.compare("symmetry-x") == 0){
                    temp_type = 8;
                }
				else if (temp.compare("parabolic") == 0) {
					temp_type = 9;
				}
				else if (temp.compare("shear") == 0) {
					temp_type = 10;
				}

                bcs.vel_type.push_back(temp_type);
                os.clear();
                os.str(std::string());

                 pbcElement = pListElement->FirstChildElement("rho_type");
                 os << pbcElement->GetText();

                temp = os.str();
                if (temp.compare("dirichlet") == 0){
                    temp_type= 1;
                }else if ( temp.compare("neumann") == 0){
                    temp_type = 2;
                }else if ( temp.compare("periodic") == 0){
                    temp_type = 3;
                }else if ( temp.compare("parabolic-N") == 0){
                    temp_type = 4;
                }else if ( temp.compare("parabolic-W") == 0){
                    temp_type= 5;
                }else if ( temp.compare("inlet") == 0){
                    temp_type = 6;
                }else if ( temp.compare("wall") == 0){
                    temp_type = 7;
                }else if ( temp.compare("symmetry-x") == 0){
                    temp_type = 8;
                }


                bcs.rho_type.push_back(temp_type);
                os.str(std::string());
                os.clear();

                //loop to next boundary condition
                pListElement = pListElement->NextSiblingElement("bc");


            }

          bcs.per_translate_x = get_xml_double("uns_boundary_conditions","periodic_translation","x",xmlDoc);
        bcs.per_translate_y = get_xml_double("uns_boundary_conditions","periodic_translation","y",xmlDoc);
        bcs.per_translate_z = get_xml_double("uns_boundary_conditions","periodic_translation","z",xmlDoc);
}

void preprocessor::parse_boundary_conditions(XMLDocument &xmlDoc, quad_bcs_plus &bcs){

    const char* parent = "boundary_conditions";



    /// WEST BCS
    bcs.w_rho = get_xml_double(parent,"west","rho",xmlDoc);
    bcs.w_u = get_xml_double(parent,"west","u",xmlDoc);
    bcs.w_v = get_xml_double(parent,"west","v",xmlDoc);
    std::string temp;

    temp = get_xml_text(parent,"west","vel_type",xmlDoc);
    if (temp.compare("dirichlet") == 0){
        bcs.w_type_vel = 1;
    }else if ( temp.compare("neumann") == 0){
        bcs.w_type_vel = 2;
    }else if ( temp.compare("periodic") == 0){
        bcs.w_type_vel = 3;
    }else if ( temp.compare("parabolic-N") == 0){
        bcs.w_type_vel = 4;
    }else if ( temp.compare("parabolic-W") == 0){
        bcs.w_type_vel = 5;
    }
	else if (temp.compare("parabolic") == 0) {
		bcs.w_type_vel = 6;
	}



    temp = get_xml_text(parent,"west","rho_type",xmlDoc);
    if (temp.compare("dirichlet") == 0){
        bcs.w_type_rho = 1;
    }else if ( temp.compare("neumann") == 0){
        bcs.w_type_rho = 2;
    }else if ( temp.compare("periodic") == 0){
        bcs.w_type_rho = 3;
    }

    /// east BCS
    bcs.e_rho = get_xml_double(parent,"east","rho",xmlDoc);
    bcs.e_u = get_xml_double(parent,"east","u",xmlDoc);
    bcs.e_v = get_xml_double(parent,"east","v",xmlDoc);


    temp = get_xml_text(parent,"east","vel_type",xmlDoc);
    if (temp.compare("dirichlet") == 0){
        bcs.e_type_vel = 1;
    }else if ( temp.compare("neumann") == 0){
        bcs.e_type_vel = 2;
    }else if ( temp.compare("periodic") == 0){
        bcs.e_type_vel = 3;
    }else if ( temp.compare("parabolic-N") == 0){
        bcs.w_type_vel = 4;
    }else if ( temp.compare("parabolic-W") == 0){
        bcs.w_type_vel = 5;
    }



    temp = get_xml_text(parent,"east","rho_type",xmlDoc);
    if (temp.compare("dirichlet") == 0){
        bcs.e_type_rho = 1;
    }else if ( temp.compare("neumann") == 0){
        bcs.e_type_rho = 2;
    }else if ( temp.compare("periodic") == 0){
        bcs.e_type_rho = 3;
    }

    /// north BCS
    bcs.n_rho = get_xml_double(parent,"north","rho",xmlDoc);
    bcs.n_u = get_xml_double(parent,"north","u",xmlDoc);
    bcs.n_v = get_xml_double(parent,"north","v",xmlDoc);


    temp = get_xml_text(parent,"north","vel_type",xmlDoc);
    if (temp.compare("dirichlet") == 0){
        bcs.n_type_vel = 1;
    }else if ( temp.compare("neumann") == 0){
        bcs.n_type_vel = 2;
    }else if ( temp.compare("periodic") == 0){
        bcs.n_type_vel = 3;
    }else if ( temp.compare("parabolic-N") == 0){
        bcs.w_type_vel = 4;
    }else if ( temp.compare("parabolic-W") == 0){
        bcs.w_type_vel = 5;
    }



    temp = get_xml_text(parent,"north","rho_type",xmlDoc);
    if (temp.compare("dirichlet") == 0){
        bcs.n_type_rho = 1;
    }else if ( temp.compare("neumann") == 0){
        bcs.n_type_rho = 2;
    }else if ( temp.compare("periodic") == 0){
        bcs.n_type_rho = 3;
    }

    /// south BCS
    bcs.s_rho = get_xml_double(parent,"south","rho",xmlDoc);
    bcs.s_u = get_xml_double(parent,"south","u",xmlDoc);
    bcs.s_v = get_xml_double(parent,"south","v",xmlDoc);


    temp = get_xml_text(parent,"south","vel_type",xmlDoc);
    if (temp.compare("dirichlet") == 0){
        bcs.s_type_vel = 1;
    }else if ( temp.compare("neumann") == 0){
        bcs.s_type_vel = 2;
    }else if ( temp.compare("periodic") == 0){
        bcs.s_type_vel = 3;
    }else if ( temp.compare("parabolic-N") == 0){
        bcs.w_type_vel = 4;
    }else if ( temp.compare("parabolic-W") == 0){
        bcs.w_type_vel = 5;
    }



    temp = get_xml_text(parent,"south","rho_type",xmlDoc);
    if (temp.compare("dirichlet") == 0){
        bcs.s_type_rho = 1;
    }else if ( temp.compare("neumann") == 0){
        bcs.s_type_rho = 2;
    }else if ( temp.compare("periodic") == 0){
        bcs.s_type_rho = 3;
    }

    // add in w velocity  afterwards for 3d
}



void preprocessor::parse_initial_conditions(XMLDocument &xmlDoc, initial_conditions &initials){

    const char* parent = "initial_conditions";
    initials.average_rho = get_xml_double(parent, "average_rho", xmlDoc);
    initials.pressure_gradient = get_xml_double(parent, "pressure_gradient", xmlDoc);
    initials.rho_gradient.x = get_xml_double(parent,"rho","gradient","x",xmlDoc);
    initials.rho_gradient.y = get_xml_double(parent,"rho","gradient","y",xmlDoc);
    initials.rho_gradient.z = get_xml_double(parent,"rho","gradient","z",xmlDoc);
    initials.rho_origin_mag.x = get_xml_double(parent,"rho","origin_magnitude","x",xmlDoc);
    initials.rho_origin_mag.y = get_xml_double(parent,"rho","origin_magnitude","y",xmlDoc);
    initials.rho_origin_mag.z = get_xml_double(parent,"rho","origin_magnitude","z",xmlDoc);
    initials.origin_loc.x = get_xml_double(parent,"origin_location","x",xmlDoc);
    initials.origin_loc.y = get_xml_double(parent,"origin_location","y",xmlDoc);
    initials.origin_loc.z = get_xml_double(parent,"origin_location","z",xmlDoc);
    initials.vel_gradient.x = get_xml_double(parent,"vel","gradient","x",xmlDoc);
    initials.vel_gradient.y = get_xml_double(parent,"vel","gradient","y",xmlDoc);
    initials.vel_gradient.z = get_xml_double(parent,"vel","gradient","z",xmlDoc);
    initials.vel_origin_mag.x = get_xml_double(parent,"vel","origin_magnitude","x",xmlDoc);
    initials.vel_origin_mag.y = get_xml_double(parent,"vel","origin_magnitude","y",xmlDoc);
    initials.vel_origin_mag.z = get_xml_double(parent,"vel","origin_magnitude","z",xmlDoc);

    initials.velocity.x = get_xml_double(parent, "u", xmlDoc);
    initials.velocity.y = get_xml_double(parent, "v", xmlDoc);
    initials.velocity.z = get_xml_double(parent, "w", xmlDoc);

	//newer approach to handle optional arguments

	XMLNode * pRoot = xmlDoc.FirstChild();
	const char * output;
	std::string temp;
	std::ostringstream os;
	int temp_type;
	if (pRoot == nullptr) return;
	XMLError eResult;

	double iOutListValue;

	XMLElement * pElement = pRoot->FirstChildElement(parent);
	if (pElement == nullptr) return;
	XMLElement * pbcElement;

	//get GPU device for multi gpu codes

	pbcElement = pElement->FirstChildElement("pressure_direction");
	if (pbcElement == nullptr) {

	}
	else {
		os << pbcElement->GetText();
		initials.pressure_direction = os.str();

		if (initials.pressure_direction == "x") {
			initials.force_x = initials.pressure_gradient;
			initials.force_y = 0.0;
			initials.force_z = 0.0;
		}
		else if (initials.pressure_direction == "y") {
			initials.force_x = 0.0;
			initials.force_y = initials.pressure_gradient;
			initials.force_z = 0.0;
		}
		else if (initials.pressure_direction == "z") {
			initials.force_x = 0.0;
			initials.force_y = 0.0;
			initials.force_z = initials.pressure_gradient;
		}
	}


}
void preprocessor::parse_geometry_variables(XMLDocument &xmlDoc, domain_geometry &geometry){

    const char* parent = "geometry";
    geometry.X = get_xml_double(parent, "x", xmlDoc);
    geometry.Y= get_xml_double(parent, "y",xmlDoc);
    geometry.Z= get_xml_double(parent, "z",xmlDoc);
    geometry.dx = get_xml_double(parent, "dx", xmlDoc);
    geometry.dy = get_xml_double(parent, "dy", xmlDoc);
    geometry.dz = get_xml_double(parent, "dz", xmlDoc);
    geometry.dt= get_xml_double(parent, "streaming_dt", xmlDoc);

	geometry.origin_x = get_xml_double(parent, "origin_x", xmlDoc);
	geometry.origin_y = get_xml_double(parent, "origin_y", xmlDoc);
	geometry.origin_z = get_xml_double(parent, "origin_z", xmlDoc);
	geometry.ref_length = get_xml_double(parent, "ref_length", xmlDoc);
	geometry.plateau_slope = get_xml_double(parent, "plateau_slope", xmlDoc);
	geometry.plateau_width = get_xml_double(parent, "plateau_width", xmlDoc);


	// newer approach to handle optional arguments

	XMLNode * pRoot = xmlDoc.FirstChild();
	const char * output;
	std::string temp;
	std::ostringstream os;
	int temp_type;
	if (pRoot == nullptr) return;
	XMLError eResult;

	double iOutListValue;

	XMLElement * pElement = pRoot->FirstChildElement(parent);
	if (pElement == nullptr) return;
	XMLElement * pbcElement;

	//get GPU device for multi gpu codes

	pbcElement = pElement->FirstChildElement("IBM_face");
	if (pbcElement == nullptr) {

	}
	else {
		os << pbcElement->GetText();
		geometry.IBM_face = os.str();

	}




    geometry.initialise();


}
void preprocessor::parse_global_variables(XMLDocument &xmlDoc, global_variables &globals){

    const char* parent = "global_variables";
	const char* target = "tolerance";
    globals.tolerance = get_xml_double(parent, target, xmlDoc);
    globals.pre_conditioned_gamma = get_xml_double(parent, "pre_condition_gamma",xmlDoc);
    globals.simulation_length = get_xml_double(parent, "simulation_length", xmlDoc);
    globals.time_marching_step = get_xml_double(parent, "time_marching_step", xmlDoc);
    globals.reynolds_number = get_xml_double(parent, "reynolds_no", xmlDoc);
    globals.max_velocity = get_xml_double(parent, "mach_no", xmlDoc);
    globals.simulation_name = get_xml_text(parent,"simulation_name",xmlDoc);
    globals.output_file_dir = get_xml_text(parent, "output_directory", xmlDoc);
    globals.max_mg_levels = get_xml_double(parent,"max_multi_grid_levels",xmlDoc);
    globals.fmg_levels =  get_xml_double(parent,"FMG_levels",xmlDoc);
    globals.fmg_tolerance = get_xml_double(parent,"FMG_tolerance",xmlDoc);
    globals.arti_disp_kappa_2 = get_xml_double(parent,"dissipation_kappa_2",xmlDoc);
    globals.arti_disp_kappa_4 = get_xml_double(parent,"dissipation_kappa_4",xmlDoc);
    globals.martinelli = get_xml_double(parent,"martinelli_exponent",xmlDoc);
    globals.testcase= get_xml_double(parent, "testcase", xmlDoc);
    globals.mesh_type= get_xml_double(parent, "mesh_type", xmlDoc);
    globals.scale= get_xml_double(parent, "scale", xmlDoc);
    globals.womersley_no= get_xml_double(parent, "womersley_no", xmlDoc);
    globals.output_step = get_xml_double(parent, "output_step", xmlDoc);
    globals.import_file = get_xml_text(parent,"import_file",xmlDoc);
    globals.import_format = get_xml_text(parent,"import_format",xmlDoc);
    globals.gradient_calc_type = get_xml_text(parent,"gradient_calc",xmlDoc);
    globals.restart_analysis = get_xml_text(parent,"restart_analysis",xmlDoc);
    globals.time_stepping = get_xml_text(parent,"time_stepping",xmlDoc);
    globals.preconditioning = get_xml_text(parent,"preconditioning",xmlDoc);
	globals.gpu_solver = get_xml_double(parent, "gpu_solver", xmlDoc);
	globals.fneq_min = get_xml_double(parent, "fneq_min", xmlDoc);

	// newer approach to handle optional arguments

	XMLNode * pRoot = xmlDoc.FirstChild();
	const char * output;
	std::string temp;
	std::ostringstream os;
	int temp_type;
	if (pRoot == nullptr) return;
	XMLError eResult;

	double iOutListValue;

	XMLElement * pElement = pRoot->FirstChildElement(parent);
	if (pElement == nullptr) return;
	XMLElement * pbcElement;

	//get GPU device for multi gpu codes

	pbcElement = pElement->FirstChildElement("gpu_device");
	if (pbcElement == nullptr) {

	}
	else {
		eResult = pbcElement->QueryDoubleText(&iOutListValue);
		globals.gpu_device = (int) iOutListValue;


	}


}

double preprocessor::get_xml_double(const char* parent, const char* child, XMLDocument &doc){

    XMLError eResult;
    XMLNode * pRoot = doc.FirstChild();

    if (pRoot == nullptr) {

        eResult = XML_ERROR_FILE_READ_ERROR;
    }

    XMLElement * parent_element = pRoot->FirstChildElement(parent);
    if (parent_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child_element = parent_element -> FirstChildElement( child );
    if (child_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;

    double output;

    eResult = child_element->QueryDoubleText(&output);

    return output;
}
double preprocessor::get_xml_double(const char* parent, const char* child, const char* child2,
                                    const char* child3, XMLDocument &doc){

    XMLError eResult;
    XMLNode * pRoot = doc.FirstChild();

    if (pRoot == nullptr) {

        eResult = XML_ERROR_FILE_READ_ERROR;
    }

    XMLElement * parent_element = pRoot->FirstChildElement(parent);
    if (parent_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child_element = parent_element -> FirstChildElement( child );
    if (child_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child2_element = child_element -> FirstChildElement( child2 );
    if (child2_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child3_element = child2_element -> FirstChildElement( child3 );
    if (child3_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    double output;

    eResult = child3_element->QueryDoubleText(&output);

    return output;
}

double preprocessor::get_xml_double(const char* parent, const char* child, const char* child2,
                                     XMLDocument &doc){

    XMLError eResult;
    XMLNode * pRoot = doc.FirstChild();

    if (pRoot == nullptr) {

        eResult = XML_ERROR_FILE_READ_ERROR;
    }

    XMLElement * parent_element = pRoot->FirstChildElement(parent);
    if (parent_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child_element = parent_element -> FirstChildElement( child );
    if (child_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child2_element = child_element -> FirstChildElement( child2 );
    if (child2_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;

    double output;

    eResult = child2_element->QueryDoubleText(&output);

    return output;
}
const char * preprocessor::get_xml_text(const char* parent, const char* child,XMLDocument &doc){

    XMLError eResult;
    XMLNode * pRoot = doc.FirstChild();

    if (pRoot == nullptr) {

        eResult = XML_ERROR_FILE_READ_ERROR;
    }


    XMLElement * parent_element = pRoot->FirstChildElement(parent);
    if (parent_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child_element = parent_element -> FirstChildElement( child );
    if (child_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;

    const char * output;

    output = child_element->GetText();

    return output;
}
const char * preprocessor::get_xml_text(const char* parent, const char* child, const char* child2
                                        , XMLDocument &doc){

    XMLError eResult;
    XMLNode * pRoot = doc.FirstChild();

    if (pRoot == nullptr) {

        eResult = XML_ERROR_FILE_READ_ERROR;
    }


    XMLElement * parent_element = pRoot->FirstChildElement(parent);
    if (parent_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child_element = parent_element -> FirstChildElement( child );
    if (child_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child2_element = child_element -> FirstChildElement( child2 );
    if (child2_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;

    const char * output;

    output = child2_element->GetText();

    return output;
}
