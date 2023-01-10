#include "Equations.h"

/**
 * @brief evaluates Schwefel’s function
 * @param vector input vector
 * @param length dimension of vector
 * @return fitness value
 */
double Equations::schwefelsfunction(double* vector, int length){

    double summation =0;

    for(int i=0; i<length; i++){
        summation += -vector[i] *sin(sqrt(fabs(vector[i])));

    }

    return 418.9829 * length -summation;

}

/**
 * @brief evaluates 1st De Jong’s function
 * @param vector input vector
 * @param length dimension of vector
 * @return fitness value
 */
double Equations::firstDeJongsfunction(double* vector, int length){

    double fitness = 0;
    for(int i=0; i< length; i++){
        fitness += vector[i] * vector[i];
    }

    return fitness;

}

/**
 * @brief evaluates Rosenbrock function
 * @param vector input vector
 * @param length dimension of vector
 * @return fitness value
 */
double Equations::rosenbrock(double* vector, int length){

    double fitness = 0;
    for(int i=0; i< length-1; i++){
        fitness += 100* pow(pow(vector[i],2) - vector[i+1] ,2) + pow(1-vector[i], 2);
    }

    return fitness;
}

/**
 * @brief evaluates Rastrigin function
 * @param vector input vector
 * @param length dimension of vector
 * @return fitness value
 */
double Equations::rastrigin(double* vector, int length){

    double summation = 0;
    for(int i=0; i< length; i++){
        summation += pow(vector[i],2) -10*cos(2*pi*vector[i]);
    }

    return 10*length*summation;
}

/**
 * @brief evaluates Griewangk function
 * @param vector input vector
 * @param length dimension of vector
 * @return fitness value
 */
double Equations::griewangk(double* vector, int length){

    double summationSection = 0;
    for(int i=0; i< length; i++){
        summationSection += pow(vector[i], 2);
    }
    summationSection /= 4000;

    double productSection = 1;
    for(int i=0; i< length; i++){
        productSection *= cos(vector[i] / sqrt(i+1));
    }

    return 1 + summationSection - productSection;
}

/**
 * @brief evaluates Sine Envelope Sine Wave function
 * @param vector input vector
 * @param length dimension of vector
 * @return fitness value
 */
double Equations::sineEnvelopeSineWave(double* vector, int length){

    double summation = 0;

    for(int i=0; i< length-1; i++){
        double numerator = pow(sin(pow(vector[i],2)+pow(vector[i+1],2) -0.5) , 2);
        double denominator = pow(1+0.001*(pow(vector[i],2)+pow(vector[i+1],2)) ,2);

        summation += 0.5 + numerator/denominator;
    }

    return -summation;

}

/**
 * @brief evaluates Stretched V Sine Wave function
 * @param vector input vector
 * @param length dimension of vector
 * @return fitness value
 */
double Equations::stretchedVSineWave(double* vector, int length){

    double fitness = 0;

    for(int i=0; i< length -1; i++){
        double firstPart = pow(pow(vector[i],2)+pow(vector[i+1],2) , 0.25);
        double secondPart = 50*pow(pow(vector[i],2)+pow(vector[i+1],2) ,0.1);

        fitness += firstPart * pow(sin(secondPart), 2)+1;
    }

    return fitness;

}

/**
 * @brief evaluates Ackley’s One function
 * @param vector input vector
 * @param length dimension of vector
 * @return fitness value
 */
double Equations::ackleysOne(double* vector, int length){

    double fitness = 0;

    for(int i=0; i< length -1; i++){
        double subsection = pow(pow(vector[i] ,2)+pow(vector[i+1] ,2) ,0.5);

        fitness += (1/exp(0.2)) * subsection + 3*(cos(2*vector[i]) + sin(2*vector[i+1]));
    }

    return fitness;

}

/**
 * @brief evaluates Ackley’s Two function
 * @param vector input vector
 * @param length dimension of vector
 * @return fitness value
 */
double Equations::ackleysTwo(double* vector, int length){

    double fitness = 0;

    for(int i=0; i< length -1; i++){


        double firstPart = cos(2*pi*vector[i])+cos(2*pi*vector[i+1]);
        double secondPart = sqrt((pow(vector[i],2)+pow(vector[i+1],2))/2);

        fitness += 20+exp(1)-20/exp(0.2*secondPart) - exp(0.5*firstPart);


    }

    return fitness;

}

/**
 * @brief evaluates Egg Holder function
 * @param vector input vector
 * @param length dimension of vector
 * @return fitness value
 */
double Equations::eggHolder(double* vector, int length){

    double fitness = 0;

    for(int i=0; i< length -1; i++){

        double firstPart = sin(sqrt(fabs(vector[i] - vector[i+1] -47)));
        double secondPart = sin(sqrt(fabs(vector[i+1]+47+vector[i]/2)));


        fitness += -vector[i]*firstPart -(vector[i+1]+47)*secondPart;

    }

    return fitness;

}

