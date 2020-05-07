#ifndef QUAD_BCS_H
#define QUAD_BCS_H


class quad_bcs
{
    public:
        quad_bcs();
        virtual ~quad_bcs();
        double w_rho, w_u,w_v,w_w;
        double s_rho, s_u, s_v, s_w;
        double e_rho, e_u, e_v, e_w;
        double n_rho, n_u, n_v, n_w;
        int w_type, s_type, e_type,n_type;
        int periodic_node;

    protected:
    private:
};

#endif // QUAD_BCS_H
