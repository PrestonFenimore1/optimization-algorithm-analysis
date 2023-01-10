//
// Created by prest on 4/6/2022.
//

#ifndef SOURCE_POPULATION_H
#define SOURCE_POPULATION_H

#include "Equations.h"
#include <chrono>
#include <random>
#include <stdexcept>

class Population{

    private:
        int rows;
        int columns;


    public:

        double ** matrix;
        double* fitness;

        Population(int experiments, int dimensionality, double lowerBound, double upperBound);
        ~Population();
        void runExperiment(int type);


};

#endif //SOURCE_POPULATION_H
