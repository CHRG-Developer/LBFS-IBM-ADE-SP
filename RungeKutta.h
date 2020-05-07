#include "Solution.h"
#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H


class RungeKutta
{
    public:
        RungeKutta();
        virtual ~RungeKutta();
//
        double alpha[4] = {1.0 , 0.5, 0.5,1.0 };
        double beta[4] = {1.0, 1.0/3.0, 1.0/3.0, 1.0/6.0} ;
        int timesteps = 1;

        //double alpha[4] = {1.0 , 0.5, 0.5,1.0 };
        //double beta[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0} ;
        //int timesteps = 4;




    protected:

    private:
};

#endif // RUNGEKUTTA_H
