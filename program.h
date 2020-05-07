#include <string>
#include "Solution.h"
#include "Mesh.h"
#include "quad_bcs_plus.h"
#include "initial_conditions.h"
#include "global_variables.h"
#include "domain_geometry.h"
#include "post_processing.h"
using namespace std;

#ifndef PROGRAM_H
#define PROGRAM_H


class program
{
    public:
        program();
        virtual ~program();
        void run (char* xml_input);
        void copyfile( char* SRC, std::string  DEST);
        void fmg_cycle(int &fmg,Solution &residual , Solution &soln,
                                      Mesh &fine_mesh, quad_bcs_plus &bcs,
                                      initial_conditions &initial_conds,
                                      global_variables globals, domain_geometry &fine_domain,
                                      Boundary_Conditions &fine_bcs , post_processing &pp);
        void output_globals (global_variables globals ,double duration,post_processing &pp);
        void output_tecplot(global_variables &globals, Mesh &Mesh, Solution &Soln,
                            Boundary_Conditions &bcs) ;
        void remove_existing_files(global_variables &globals);
        void remove_directory_files(std::string output_location);
    protected:
    private:
};

#endif // PROGRAM_H
