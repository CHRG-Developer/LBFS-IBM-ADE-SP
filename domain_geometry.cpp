#include "domain_geometry.h"
#include <math.h>

domain_geometry::domain_geometry()
{
    //ctor
}

domain_geometry::~domain_geometry()
{
    //dtor
}

void domain_geometry::initialise(){

    cs = 1/sqrt(3);


}

void domain_geometry::scale_geometries(double scale){
    X = X*scale ;
    Y = Y*scale ;

}
