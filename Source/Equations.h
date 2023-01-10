#ifndef SOURCE_EQUATIONS_H
#define SOURCE_EQUATIONS_H

#include <cmath>
#define pi 4*atan(1)

class Equations{

    public:
        double schwefelsfunction(double* vector, int length);
        double firstDeJongsfunction(double* vector, int length);
        double rosenbrock(double* vector, int length);
        double rastrigin(double* vector, int length);
        double griewangk(double* vector, int length);
        double sineEnvelopeSineWave(double* vector, int length);
        double stretchedVSineWave(double* vector, int length);
        double ackleysOne(double* vector, int length);
        double ackleysTwo(double* vector, int length);
        double eggHolder(double* vector, int length);

};
typedef double (Equations::*EquationsMemFn)(double* vector, int length);

#endif //SOURCE_EQUATIONS_H
