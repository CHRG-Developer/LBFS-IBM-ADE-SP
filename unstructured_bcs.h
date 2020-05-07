#ifndef UNSTRUCTURED_BCS_H
#define UNSTRUCTURED_BCS_H

#include <vector>
#include <string>
class unstructured_bcs
{
    public:
        unstructured_bcs();
        virtual ~unstructured_bcs();

        std::vector<std::string> name;
         std::vector<int> rho_type, vel_type;
        std::vector<double> rho,u,v,w;

        double per_translate_x,per_translate_y, per_translate_z;

    protected:

    private:
};

#endif // UNSTRUCTURED_BCS_H
