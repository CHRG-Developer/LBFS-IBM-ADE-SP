#ifndef PREPROCESSOR_H
#define PREPROCESSOR_H
#include "global_variables.h"
#include "tinyxml2.h"
#include "domain_geometry.h"
#include "initial_conditions.h"
#include "quad_bcs_plus.h"
#include "unstructured_bcs.h"
#include <vector>
#include "lagrangian_object.h"

using namespace tinyxml2;

class preprocessor
{
    public:
        preprocessor();
        virtual ~preprocessor();
        void initialise_program_variables(char* xml_input, global_variables &globals,
                                          domain_geometry &geometry,initial_conditions &initial_conds,
                                           quad_bcs_plus &bcs,unstructured_bcs &u_bcs , std::vector <lagrangian_object> &object_vec);

    protected:
    private:
        void parse_global_variables(XMLDocument &xmlDoc, global_variables &globals);
        void parse_geometry_variables(XMLDocument &xmlDoc, domain_geometry &geometry);
        void parse_initial_conditions(XMLDocument &xmlDoc, initial_conditions &initials);
        void parse_boundary_conditions(XMLDocument &xmlDoc, quad_bcs_plus &bcs);

        double get_xml_double(const char* parent, const char* child, XMLDocument &doc);
        double get_xml_double(const char* parent, const char* child, const char* child2,
                                    const char* child3, XMLDocument &doc);
        double get_xml_double(const char* parent, const char* child, const char* child2,
                                     XMLDocument &doc);
        const char * get_xml_text(const char* parent, const char* child, XMLDocument &doc);
        const char * get_xml_text(const char* parent, const char* child, const char* child2
                                        , XMLDocument &doc);
        void mach_number_factor( global_variables &globals,quad_bcs_plus &bcs,
        initial_conditions &initials,domain_geometry &geometry );

        void parse_uns_boundary_conditions(XMLDocument &xmlDoc, unstructured_bcs &bcs
                ,global_variables &globals);
		void parse_lagrangian_objects(XMLDocument &xmlDoc, std::vector <lagrangian_object> &object_vec, global_variables &globals);
};

#endif // PREPROCESSOR_H
