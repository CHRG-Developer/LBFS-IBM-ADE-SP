#include "external_forces.h"
#include <stdlib.h>
using namespace std;
external_forces::external_forces(int _total_nodes)
{
    //ctor
    total_nodes = _total_nodes;

    force = new double [total_nodes ];
    if (force==NULL) exit (1);
}

external_forces::~external_forces()
{
    //dtor
    delete [] force;
    force = NULL;
}


void external_forces::set_uniform_force(double magnitude){

    for( int i =0; i< total_nodes; i++){
        force[i] = magnitude;

    }

}
