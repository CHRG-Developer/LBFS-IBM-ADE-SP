#ifndef QUAD_BCS_PLUS_H
#define QUAD_BCS_PLUS_H


class quad_bcs_plus
{
    public:
        quad_bcs_plus();
        virtual ~quad_bcs_plus();
        double w_rho, w_u,w_v,w_w;
        double s_rho, s_u, s_v, s_w;
        double e_rho, e_u, e_v, e_w;
        double n_rho, n_u, n_v, n_w;
        int w_type_rho, s_type_rho, e_type_rho,n_type_rho;
        int w_type_vel, s_type_vel, e_type_vel,n_type_vel;
        int periodic_node;

    protected:
    private:
};

#endif // QUAD_BCS_H
