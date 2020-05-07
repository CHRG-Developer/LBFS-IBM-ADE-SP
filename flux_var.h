#ifndef FLUX_VAR_H
#define FLUX_VAR_H


class flux_var
{
    public:
        flux_var();
        virtual ~flux_var();

        double P;
        double momentum_x;
        double momentum_y;
        double momentum_z;

        void div_volume( double vol, double small_num);
        void zero(){P = 0.0; momentum_x = 0.0; momentum_y =0.0; momentum_z =0.0;};
    protected:

    private:
};

#endif // FLUX_VAR_H
