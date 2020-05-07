#ifndef DOMAIN_GEOMETRY_H
#define DOMAIN_GEOMETRY_H
#include <string>

class domain_geometry
{
    public:
        domain_geometry();
        virtual ~domain_geometry();
        double X,Y,Z,dx,dy,dz;
        double dt ;  // streaming time step
        double cs; // speed of sound in medium
		double ref_length, origin_x, origin_y, origin_z;
        void initialise();
        void scale_geometries(double scale);
		double plateau_width, plateau_slope;
		std::string IBM_face;

    protected:
    private:
};

#endif // DOMAIN_GEOMETRY_H
