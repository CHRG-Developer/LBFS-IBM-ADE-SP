#include "vector_var.h"
#include <math.h>       /* pow ,cosh */

vector_var::vector_var()
{
    //ctor
}
vector_var::vector_var(double _x, double _y, double _z)
{
    //ctor
    x = _x;
    y = _y;
    z = _z;
}
//vector_var::~vector_var()
//{
//    //dtor
//}


void vector_var::zero()
{
    //ctor
    x = 0.0;
    y = 0.0;
    z = 0.0;
}

double  vector_var::Dot_Product(vector_var b){

    return x * b.x + y *b.y + z * b.z;

}

void  vector_var::cross_product(vector_var b){
    double x1,y1,z1;

    x1 = y * b.z - z*b.y;
    y1 = z* b.x - x*b.z;
    z1 = x*b.y - y*b.x;

    x = x1;
    y = y1;
    z = z1;
}
void  vector_var::set_equal(vector_var b){

    x = b.x;
    y = b.y;
    z = b.z;
}
void  vector_var::add(vector_var b){

    x = x + b.x;
    y = y + b.y;
    z = z + b.z;
}
void  vector_var::average(vector_var b,vector_var c,vector_var d){

    x =  (b.x + c.x + d.x)/3;
    y = (b.y + c.y + d.y)/3;
    z = (b.z + c.z + d.z)/3;
}

void  vector_var::subtract(vector_var b){

    x = x - b.x;
    y = y - b.y;
    z = z - b.z;
}


vector_var vector_var::line_magnitude(vector_var intercept, vector_var slope, vector_var displacement){

    vector_var temp;

    temp.x = intercept.x + slope.x * displacement.x;
    temp.y = intercept.y + slope.y *displacement.y;
    temp.z = intercept.z + slope.z * displacement.z;

    return temp;


}
double vector_var::Magnitude(){

    return sqrt(pow(x,2) + pow(y,2) + pow(z,2) );
}

double  vector_var::Angle_Between_Vectors(vector_var b){
    double num, denom;
    num = Dot_Product(b);
    denom = Magnitude();
    denom = denom* b.Magnitude();
    num = num/denom;
    num = acos (num);

    return acos (Dot_Product(b)/(Magnitude()* b.Magnitude())) ;

}

void vector_var::factor(double a){

    x = x/a;
    y= y/a;
    z= z/a;
}
void vector_var::round(int precision){



    double f = pow(10,precision);
    double temp;


    temp = floor(fabs(x *f))/f;
    if (x < 0.0){
        x = -temp;
    }else{
        x = temp;
    }
    temp = floor(fabs(y *f))/f;
    if (y < 0.0){
        y = -temp;
    }else{
        y = temp;
    }
    temp = floor(fabs(z *f))/f;
    if (z < 0.0){
        z = -temp;
    }else{
        z = temp;
    }

}
void vector_var::Get_Gradient(double y1, double y2, vector_var x1, vector_var x2 ){

        x = y2 -y1 ;
        y = y2 -y1 ;

        //check if x1 and x2 are same co-ordinate -> zero gradient
        if( fabs(x2.x - x1.x ) < pow(10,-5) ){
            x = 0;
        }else{
            x =  x / (x2.x - x1.x );
        }
        if( fabs(x2.y -x1.y) < pow(10,-5) ){
            y = 0;
        }else{
            y = y / (x2.y -x1.y);
        }
        z = 0; //update for 3d


}

void vector_var::relative_vectors(vector_var origin, vector_var ref1, vector_var ref2, double const2){

    x = origin.x -ref1.x - ref2.x * const2;
    y = origin.y -ref1.y - ref2.y * const2;
    z = 0; // update in 3d


}
