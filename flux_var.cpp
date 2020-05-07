#include "flux_var.h"

flux_var::flux_var()
{
    //ctor
}

flux_var::~flux_var()
{
    //dtor
}

void flux_var::div_volume(double vol, double small_num){
    
    if( vol > small_num)
    P = P/vol;
    momentum_x = momentum_x/vol;
    momentum_y = momentum_y/vol;
    momentum_z = momentum_z/vol;
    
    
}