//
// Created by prest on 4/19/2022.
//

#ifndef SOURCE_ALGORITHMS_H
#define SOURCE_ALGORITHMS_H

#include "Population.h"
#include "Equations.h"
#include <functional>
#include <stdexcept>

class Algorithms {

    int equationNum;
    int dimension;
    double lowerBound;
    double upperBound;
    Equations eq;
    EquationsMemFn equation;

    private:

        void addVectors(double* result, double* arr1, double* arr2, int length);
        void subtractVectors(double* result, double* arr1, double* arr2, int length);
        void multiplyByScalar(double* result, double scalar, int length);
        void mutuallyExclusiveParents(int parents[], int numToMake, int disallowedValue, std::mt19937* g1, std::uniform_int_distribution<>* d);
        void forceInBounds(double* vector);

        void deMutationStrategy1or6(Population* pop, double* popFitness, double** noisyVectors, int populationSize,
            double scalingFactor, double lambda, std::mt19937* g1);
        void deMutationStrategy2or7(Population* pop, double* popFitness, double** noisyVectors, int populationSize,
            double scalingFactor, double lambda, std::mt19937* g1);
        void deMutationStrategy3or8(Population* pop, double* popFitness, double** noisyVectors, int populationSize,
            double scalingFactor, double lambda, std::mt19937* g1);
        void deMutationStrategy4or9(Population* pop, double* popFitness, double** noisyVectors, int populationSize,
            double scalingFactor, double lambda, std::mt19937* g1);
        void deMutationStrategy5or10(Population* pop, double* popFitness, double** noisyVectors, int populationSize,
            double scalingFactor, double lambda, std::mt19937* g1);

        void exponentialCrossover(Population* pop, double* popFitness, double** noisyVectors, int populationSize, double crossoverFactor, std::mt19937* g1);
        void binomialCrossover(Population* pop, double* popFitness, double** noisyVectors, int populationSize, double crossoverFactor, std::mt19937* g1);



    public:

        double* bestFit;

        Algorithms(int equationNumber, int dimension, double lowerBound, double upperBound);
        ~Algorithms();

        double randomWalk(int iterations);
        double localSearch();
        double repeatedLocalSearch(int iterations);
        double differentialEvolution(int populationSize, double crossoverFactor, double scalingFactor, double lambda, int generations, int strategy);
        double particleSwarmOptimization(int particles, double c1, double c2, int generations);

};
typedef void (Algorithms::*mutationFunction)(Population* pop, double* popFitness, double** noisyVectors,
        int populationSize, double scalingFactor, double lambda, std::mt19937* g1);
typedef void (Algorithms::*crossoverFunction)(Population* pop, double* popFitness, double** noisyVectors,
        int populationSize, double crossoverFactor, std::mt19937* g1);

#endif //SOURCE_ALGORITHMS_H
