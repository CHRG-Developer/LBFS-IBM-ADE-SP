#ifndef EXTERNAL_FORCES_H
#define EXTERNAL_FORCES_H


class external_forces
{
    public:
        external_forces(int _total_nodes);
        virtual ~external_forces();
        void set_uniform_force(double magnitude);
        double get_force( int i) {return force[i];};
    protected:
    private:
        double *force;
        int total_nodes;
};

#endif // EXTERNAL_FORCES_H
