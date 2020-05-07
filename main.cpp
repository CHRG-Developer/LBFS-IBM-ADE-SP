#include <iostream>
#include "Mesh.h"
#include "vector_var.h"
#include <stdio.h>      /* printf */
#include <iostream>
#include <math.h>
#include "Boundary_Conditions.h"
#include "Solution.h"
#include "Solver.h"
#include "test_cases.h"
#include "program.h"




using namespace std;


int main(int argc, char *argv[])
{
    // Xml input file from command line argument
    char *xml_input;
    xml_input = argv[1] ;



    program LBFS;

    LBFS.run(xml_input);




    //test_cases lid_driven_cav;
    //lid_driven_cav.lid_driven_cavity_N();
   //  couette_flow.west_to_east_couette_flow();
 //   couette_flow.east_to_west_couette_flow();
    // couette_flow.north_to_south_couette_flow();
   // couette_flow.south_to_north_couette_flow();

    return 0;
};

