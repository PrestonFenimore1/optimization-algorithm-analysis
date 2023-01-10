//
// Created by prest on 4/19/2022.
//

#include "Algorithms.h"

/**
 * @param equationNumber which equation will be considered for class methods
 * @param dimension the dimension of the vectors to be used
 * @param lowerBound the lower bound of the solution space
 * @param upperBound the upper bound of the solution space
 * @throws invalid_argument when equationNumber > 10 or <1 or dimension < 1
 */
Algorithms::Algorithms(int equationNumber, int dimension, double lowerBound, double upperBound){

    if(equationNumber > 10 || equationNumber < 1){
        throw std::invalid_argument( "equationNumber must be 1-10" );
    }

    if(dimension < 2){
        throw std::invalid_argument( "dimension must be greater than or equal to 2" );
    }

    Algorithms::equationNum = equationNumber;
    Algorithms::lowerBound = lowerBound;
    Algorithms::upperBound = upperBound;
    Algorithms::dimension = dimension;
    Algorithms::bestFit = new double[dimension];

    switch(equationNum){
        case 1:
            Algorithms::equation = &Equations::schwefelsfunction;
            break;
        case 2:
            Algorithms::equation = &Equations::firstDeJongsfunction;
            break;
        case 3:
            Algorithms::equation = &Equations::rosenbrock;
            break;
        case 4:
            Algorithms::equation = &Equations::rastrigin;
            break;
        case 5:
            Algorithms::equation = &Equations::griewangk;
            break;
        case 6:
            Algorithms::equation = &Equations::sineEnvelopeSineWave;
            break;
        case 7:
            Algorithms::equation = &Equations::stretchedVSineWave;
            break;
        case 8:
            Algorithms::equation = &Equations::ackleysOne;
            break;
        case 9:
            Algorithms::equation = &Equations::ackleysTwo;
            break;
        case 10:
            Algorithms::equation = &Equations::eggHolder;
            break;
    }

}

Algorithms::~Algorithms(){
    delete[] Algorithms::bestFit;
}

/**
 * performs the random walk algorithm
 * stores best vector in the bestVector class variable
 * @param iterations how many vectors will by analyzed
 * @return best fitness found
 * @throws invalid_argument when iterations < 1
 */
double Algorithms::randomWalk(int iterations){

    if(iterations < 1){
        throw std::invalid_argument( "received value less than 1" );
    }

    //generates a set of solution vectors
    Population pop(iterations, Algorithms::dimension,Algorithms::lowerBound, Algorithms::upperBound);

    double* localBestFit = pop.matrix[0];

    double bestFitness = std::invoke(Algorithms::equation, Algorithms::eq, localBestFit, Algorithms::dimension);

    //checks fitness of every generated solution vector
    for(int i=1; i<iterations; i++){

        double newFitness = std::invoke(Algorithms::equation, Algorithms::eq, pop.matrix[i], Algorithms::dimension);

        if(newFitness < bestFitness){
            bestFitness = newFitness;
            localBestFit = pop.matrix[i];
        }
    }

    for(int i=0; i<Algorithms::dimension; i++){
        Algorithms::bestFit[i] = localBestFit[i];
    }

    return bestFitness;

}

/**
 * performs the local search algorithm
 * stores best vector in the bestVector class variable
 * @return best fitness found
 */
double Algorithms::localSearch(){

    //create random initial vector
    Population pop(1, Algorithms::dimension, Algorithms::lowerBound, Algorithms::upperBound);

    //copy random initial vector to class variable
    for(int i=0; i<Algorithms::dimension; i++){
        Algorithms::bestFit[i] = pop.matrix[0][i];
    }

    //get fitness of initial vector
    double bestValue = std::invoke(Algorithms::equation, Algorithms::eq, bestFit, Algorithms::dimension);
    //allocate memory for result vector from gradient decent
    double* testVector = new double[Algorithms::dimension];
    double alpha = 0.3;

    bool foundBetterSolution = true;


    while(foundBetterSolution){

        for(int i=0; i<Algorithms::dimension; i++){

            Algorithms::bestFit[i] += alpha;
            double testFitness = std::invoke(Algorithms::equation, Algorithms::eq, Algorithms::bestFit, Algorithms::dimension);
            Algorithms::bestFit[i] -= alpha;

            testVector[i] = Algorithms::bestFit[i] - alpha*(testFitness - bestValue);

            //assures testVector is within solution space
            if(testVector[i] < Algorithms::lowerBound){
                testVector[i] = Algorithms::lowerBound;
            }
            else if(testVector[i] > Algorithms::upperBound){
                testVector[i] = Algorithms::upperBound;
            }

        }

        //calculates fitness of obtained vector
        double newFitness = std::invoke(Algorithms::equation, Algorithms::eq, testVector, Algorithms::dimension);

        //checks if new vector is better than old one
        //if not, loop is terminated, else loop is restarted with new vector as initial values
        if(bestValue <= newFitness){
            foundBetterSolution = false;
        }
        else{
            for(int i=0; i<Algorithms::dimension; i++){
                Algorithms::bestFit[i] = testVector[i];
                pop.matrix[0][i] = testVector[i];
            }
            bestValue = newFitness;
        }

    }

    delete[] testVector;

    return bestValue;
}

/**
 * performs the repeated local search algorithm
 * stores best vector in the bestVector class variable
 * @param iterations how many local searches will be performed
 * @return best fitness found
 * @throws invalid_argument when iterations < 1
 */
double Algorithms::repeatedLocalSearch(int iterations){

    if(iterations < 1){
        throw std::invalid_argument( "received value less than 1" );
    }

    //sets base case for algorithm
    double bestValue = localSearch();

    //sets generated best fitness vector as local variable to serve as best found vector
    double* localBestFit = new double[Algorithms::dimension];
    for(int i=0; i<Algorithms::dimension; i++){
        localBestFit[i] = Algorithms::bestFit[i];
    }

    for(int i=1; i<iterations; i++){

        double testValue = localSearch();

        //checks of better vector was found
        if(testValue < bestValue){

            for(int j=0; j<Algorithms::dimension; j++){
                localBestFit[j] = Algorithms::bestFit[j];
            }
            bestValue = testValue;

        }

    }

    //set class variable to match local variable
    for(int i=0; i<Algorithms::dimension; i++){
        Algorithms::bestFit[i] = localBestFit[i];
    }

    delete[] localBestFit;

    return bestValue;

}


/**
 * adds two vectors and saves result to first vector given as parameter
 * @param result pointer to resultant vector
 * @param arr1
 * @param arr2
 * @param length length of vectors
 */
void Algorithms::addVectors(double* result, double* arr1, double* arr2, int length){
    for(int i=0; i<length; i++){
        result[i] = arr1[i]+arr2[i];
    }
}

/**
 * subtracts arr2 from arr1 and saves result to first vector given as parameter
 * @param result pointer to resultant vector
 * @param arr1
 * @param arr2
 * @param length length of vectors
 */
void Algorithms::subtractVectors(double* result, double* arr1, double* arr2, int length){
    for(int i=0; i<length; i++){
        result[i] = arr1[i]-arr2[i];
    }
}

/**
 * multiplies given vector by scalar
 * @param result pointer to resultant vector
 * @param arr
 * @param scalar
 * @param length length of vector
 */
void Algorithms::multiplyByScalar(double* result, double scalar, int length){
    for(int i=0; i<length; i++){
        result[i] = scalar*result[i];
    }
}

/**
 * helper function for mutation strategies
 * creates mutually exclusive parents for mutation
 * @param parents
 * @param numToMake
 * @param disallowedValue the current population vector for which a noisy vector is being created for
 * @param g1
 * @param d
 */
void Algorithms::mutuallyExclusiveParents(int parents[], int numToMake, int disallowedValue, std::mt19937* g1, std::uniform_int_distribution<>* d){

    for(int j=0; j<numToMake;){
        int trialParent =  (*d)(*g1);

        bool canBeAssigned = true;

        if(trialParent == disallowedValue){
            canBeAssigned = false;
        }
        else{
            for(int k=0; k<j; k++){

                if(trialParent == parents[k]){
                    canBeAssigned = false;
                }
            }
        }

        if(canBeAssigned){
            parents[j] = trialParent;
            j++;
        }
    }

}

/**
 * forces vector elements to be within solution space
 * @param noisyVector
 */
void Algorithms::forceInBounds(double* vector){

    for(int i=0; i<Algorithms::dimension; i++){

        double value = vector[i];

        if(value < Algorithms::lowerBound){
            vector[i]= Algorithms::lowerBound;
        }
        else if(value > Algorithms::upperBound){
            vector[i] = Algorithms::upperBound;
        }
    }
}

/**
 * mutation algorithm used for Differential Evolution strategies 1 and 6
 * @param pop
 * @param popFitness vector containing the fitness values for each vector in pop
 * @param noisyVectors 2d array containing arrays for each noisy vector generated for each vector in population
 * @param populationSize
 * @param scalingFactor
 * @param lambda unused, but needed for compatibility with std::invoke
 * @param g1 random generator
 */
void Algorithms::deMutationStrategy1or6(Population* pop, double* popFitness, double** noisyVectors,
    int populationSize, double scalingFactor, double lambda, std::mt19937* g1){

    //find best vector in population
    double bestFitness = popFitness[0];
    int bestVector = 0;
    for(int i=1; i< populationSize; i++){
        double newBestFitness = popFitness[i];

        if(newBestFitness < bestFitness){
            bestFitness = newBestFitness;
            bestVector = i;
        }
    }

    std::uniform_int_distribution<> d(0,populationSize-1);
    int parents[2] = {-1,-1};

    for(int i=0; i<populationSize; i++){

        //assign mutually exclusive parents
        Algorithms::mutuallyExclusiveParents(parents, 2, i, g1, &d);

        //calculate noisy vector
        Algorithms::subtractVectors(noisyVectors[i], pop->matrix[parents[0]], pop->matrix[parents[1]], Algorithms::dimension);
        Algorithms::multiplyByScalar(noisyVectors[i], scalingFactor, Algorithms::dimension);
        Algorithms::addVectors(noisyVectors[i], pop->matrix[bestVector], noisyVectors[i], Algorithms::dimension);
        forceInBounds(noisyVectors[i]);

    }

}

/**
 * mutation algorithm used for Differential Evolution strategies 2 and 7
 * @param pop
 * @param popFitness unused, but needed for compatibility with std::invoke
 * @param noisyVectors 2d array containing arrays for each noisy vector generated for each vector in population
 * @param populationSize
 * @param scalingFactor
 * @param lambda unused, but needed for compatibility with std::invoke
 * @param g1 random generator
 */
void Algorithms::deMutationStrategy2or7(Population* pop, double* popFitness, double** noisyVectors,
    int populationSize, double scalingFactor, double lambda, std::mt19937* g1 ){

    std::uniform_int_distribution<> d(0,populationSize-1);
    int parents[3] = {-1,-1, -1};

    for(int i=0; i<populationSize; i++){

        //assign mutually exclusive parents
        Algorithms::mutuallyExclusiveParents(parents, 3, i, g1, &d);

        //calculate noisy vector
        Algorithms::subtractVectors(noisyVectors[i], pop->matrix[parents[1]], pop->matrix[parents[2]], Algorithms::dimension);
        Algorithms::multiplyByScalar(noisyVectors[i], scalingFactor, Algorithms::dimension);
        Algorithms::addVectors(noisyVectors[i], pop->matrix[parents[0]], noisyVectors[i], Algorithms::dimension);
        forceInBounds(noisyVectors[i]);

    }

}

/**
 * mutation algorithm used for Differential Evolution strategies 3 and 8
 * @param pop
 * @param popFitness vector containing the fitness values for each vector in pop
 * @param noisyVectors 2d array containing arrays for each noisy vector generated for each vector in population
 * @param populationSize
 * @param scalingFactor
 * @param lambda second scaling factor
 * @param g1 random generator
 */
void Algorithms::deMutationStrategy3or8(Population* pop, double* popFitness, double** noisyVectors,
    int populationSize, double scalingFactor, double lambda, std::mt19937* g1){

    //find best vector in population
    double bestFitness = popFitness[0];
    int bestVector = 0;
    for(int i=1; i< populationSize; i++){
        double newBestFitness = popFitness[i];

        if(newBestFitness < bestFitness){
            bestFitness = newBestFitness;
            bestVector = i;
        }
    }

    std::uniform_int_distribution<> d(0,populationSize-1);
    int parents[2] = {-1,-1};

    double* tempVector = new double[Algorithms::dimension];
    for(int i=0; i<populationSize; i++){

        //assign mutually exclusive parents
        Algorithms::mutuallyExclusiveParents(parents, 2, i, g1, &d);

        //calculate noisy vector
        Algorithms::subtractVectors(noisyVectors[i], pop->matrix[parents[0]], pop->matrix[parents[1]], Algorithms::dimension);
        Algorithms::multiplyByScalar(noisyVectors[i], scalingFactor, Algorithms::dimension);

        Algorithms::subtractVectors(tempVector, pop->matrix[bestVector], pop->matrix[i], Algorithms::dimension);
        Algorithms::multiplyByScalar(tempVector, lambda, Algorithms::dimension);

        Algorithms::addVectors(noisyVectors[i], tempVector, noisyVectors[i], Algorithms::dimension);
        Algorithms::addVectors(noisyVectors[i], pop->matrix[i], noisyVectors[i], Algorithms::dimension);
        forceInBounds(noisyVectors[i]);

    }

    delete[] tempVector;
}

/**
 * mutation algorithm used for Differential Evolution strategies 4 and 9
 * @param pop
 * @param popFitness vector containing the fitness values for each vector in pop
 * @param noisyVectors 2d array containing arrays for each noisy vector generated for each vector in population
 * @param populationSize
 * @param scalingFactor
 * @param lambda unused, but needed for compatibility with std::invoke
 * @param g1 random generator
 */
void Algorithms::deMutationStrategy4or9(Population* pop, double* popFitness, double** noisyVectors,
    int populationSize, double scalingFactor, double lambda, std::mt19937* g1){

    //find best vector in population
    double bestFitness = popFitness[0];
    int bestVector = 0;
    for(int i=1; i< populationSize; i++){
        double newBestFitness = popFitness[i];

        if(newBestFitness < bestFitness){
            bestFitness = newBestFitness;
            bestVector = i;
        }
    }

    std::uniform_int_distribution<> d(0,populationSize-1);
    int parents[4] = {-1,-1, -1, -1};

    for(int i=0; i<populationSize; i++){

        //assign mutually exclusive parents
        Algorithms::mutuallyExclusiveParents(parents, 4, i, g1, &d);

        //calculate noisy vector
        Algorithms::addVectors(noisyVectors[i], pop->matrix[parents[0]], pop->matrix[parents[1]], Algorithms::dimension);
        Algorithms::subtractVectors(noisyVectors[i], noisyVectors[i], pop->matrix[parents[2]], Algorithms::dimension);
        Algorithms::subtractVectors(noisyVectors[i], noisyVectors[i], pop->matrix[parents[3]], Algorithms::dimension);
        Algorithms::multiplyByScalar(noisyVectors[i], scalingFactor, Algorithms::dimension);
        Algorithms::addVectors(noisyVectors[i], noisyVectors[i], pop->matrix[bestVector], Algorithms::dimension);
        forceInBounds(noisyVectors[i]);

    }

}

/**
 * mutation algorithm used for Differential Evolution strategies 5 and 10
 * @param pop
 * @param popFitness unused, but needed for compatibility with std::invoke
 * @param noisyVectors 2d array containing arrays for each noisy vector generated for each vector in population
 * @param populationSize
 * @param scalingFactor
 * @param lambda unused, but needed for compatibility with std::invoke
 * @param g1 random generator
 */
void Algorithms::deMutationStrategy5or10(Population* pop, double* popFitness, double** noisyVectors,
    int populationSize, double scalingFactor, double lambda, std::mt19937* g1){

    std::uniform_int_distribution<> d(0,populationSize-1);
    int parents[5] = {-1,-1, -1, -1, -1};

    for(int i=0; i<populationSize; i++){

        //assign mutually exclusive parents
        Algorithms::mutuallyExclusiveParents(parents, 5, i, g1, &d);


        //calculate noisy vector
        Algorithms::addVectors(noisyVectors[i], pop->matrix[parents[0]], pop->matrix[parents[1]], Algorithms::dimension);
        Algorithms::subtractVectors(noisyVectors[i], noisyVectors[i], pop->matrix[parents[2]], Algorithms::dimension);
        Algorithms::subtractVectors(noisyVectors[i], noisyVectors[i], pop->matrix[parents[3]], Algorithms::dimension);
        Algorithms::multiplyByScalar(noisyVectors[i], scalingFactor, Algorithms::dimension);
        Algorithms::addVectors(noisyVectors[i], noisyVectors[i], pop->matrix[parents[4]], Algorithms::dimension);
        forceInBounds(noisyVectors[i]);

    }


}

/**
 * helper function for differentialEvolution that performs exponential crossovers
 * @param pop
 * @param popFitness vector containing the fitness values for each vector in pop
 * @param noisyVectors 2d array containing arrays for each noisy vector generated for each vector in population
 * @param populationSize
 * @param crossoverFactor
 * @param g1 random generator
 */
void Algorithms::exponentialCrossover(Population* pop, double* popFitness, double** noisyVectors, int populationSize, double crossoverFactor, std::mt19937* g1){

    //not exponentialCrossover

    std::uniform_real_distribution<> d(0,1);
    std::uniform_int_distribution<> startGenerator(0, Algorithms::dimension-1);

    for(int i=0; i<populationSize; i++){

        int start = startGenerator(*g1);

        /*noisyVectors memory location is used for the trial vector
         * so every element until the starting point must be replaced
         * with the population vector elements*/
        for(int j=0; j<start; j++){
            noisyVectors[i][j] = pop->matrix[i][j];
        }

        //skip elements that were to be crossed over
        while( d(*g1) < crossoverFactor && start < Algorithms::dimension){
            start++;
        }

        //fill out rest of elements with pop elements
        for(;start<Algorithms::dimension; start++){
            noisyVectors[i][start] = pop->matrix[i][start];
        }

        double trialFitness = std::invoke(Algorithms::equation, Algorithms::eq, noisyVectors[i], Algorithms::dimension);

        //copy trial vector to population if trial vector has better fitness
        if(trialFitness < popFitness[i]){
            popFitness[i] = trialFitness;

            for(int j=0; j<Algorithms::dimension; j++){
                pop->matrix[i][j] = noisyVectors[i][j];
            }

        }

    }

}

/**
 * helper function for differentialEvolution that performs binomial crossovers
 * @param pop
 * @param popFitness vector containing the fitness values for each vector in pop
 * @param noisyVectors 2d array containing arrays for each noisy vector generated for each vector in population
 * @param populationSize
 * @param crossoverFactor
 * @param g1 random generator
 */
void Algorithms::binomialCrossover(Population* pop, double* popFitness, double** noisyVectors, int populationSize, double crossoverFactor, std::mt19937* g1){

    std::uniform_real_distribution<> d(0,1);
    std::uniform_int_distribution<> startGenerator(0, Algorithms::dimension-1);

    for(int i=0; i<populationSize; i++){


        int start = startGenerator(*g1);

        /*noisyVectors memory location is used for the trial vector
         * so every element until the starting point must be replaced
         * with the population vector elements*/
        for(int j=0; j<start; j++){
            noisyVectors[i][j] = pop->matrix[i][j];
        }

        //calculate trial vector
        for(int j=start; j<Algorithms::dimension; j++){

            double rand = d(*g1);

            //inverted because of source destination
            if(rand >= crossoverFactor){
                noisyVectors[i][j] = pop->matrix[i][j];
            }
        }

        double trialFitness = std::invoke(Algorithms::equation, Algorithms::eq, noisyVectors[i], Algorithms::dimension);

        //copy trial vector to population if trial vector has better fitness
        if(trialFitness < popFitness[i]){
            popFitness[i] = trialFitness;

            for(int j=0; j<Algorithms::dimension; j++){
                pop->matrix[i][j] = noisyVectors[i][j];
            }

        }

    }

}

/**
 * performs differential evolution
 * @param populationSize
 * @param crossoverFactor
 * @param scalingFactor
 * @param lambda a second scaling factor used in some mutation strategy
 * @param generations maximum generations
 * @param strategy mutation strategy
 * @return best fitness found
 */
double Algorithms::differentialEvolution(int populationSize, double crossoverFactor, double scalingFactor, double lambda, int generations, int strategy){

    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 g1(seed1);//mersenne_twister_engine
    Population pop(populationSize, Algorithms::dimension, Algorithms::lowerBound, Algorithms::upperBound);
    mutationFunction mutationStrategy;
    crossoverFunction crossoverStrategy;

    switch(strategy){
        case 1:
            mutationStrategy = &Algorithms::deMutationStrategy1or6;
            crossoverStrategy = &Algorithms::exponentialCrossover;
            break;
        case 2:
            mutationStrategy = &Algorithms::deMutationStrategy2or7;
            crossoverStrategy = &Algorithms::exponentialCrossover;
            break;
        case 3:
            mutationStrategy = &Algorithms::deMutationStrategy3or8;
            crossoverStrategy = &Algorithms::exponentialCrossover;
            break;
        case 4:
            mutationStrategy = &Algorithms::deMutationStrategy4or9;
            crossoverStrategy = &Algorithms::exponentialCrossover;
            break;
        case 5:
            mutationStrategy = &Algorithms::deMutationStrategy5or10;
            crossoverStrategy = &Algorithms::exponentialCrossover;
            break;
        case 6:
            mutationStrategy = &Algorithms::deMutationStrategy1or6;
            crossoverStrategy = &Algorithms::binomialCrossover;
            break;
        case 7:
            mutationStrategy = &Algorithms::deMutationStrategy2or7;
            crossoverStrategy = &Algorithms::binomialCrossover;
            break;
        case 8:
            mutationStrategy = &Algorithms::deMutationStrategy3or8;
            crossoverStrategy = &Algorithms::binomialCrossover;
            break;
        case 9:
            mutationStrategy = &Algorithms::deMutationStrategy4or9;
            crossoverStrategy = &Algorithms::binomialCrossover;
            break;
        case 10:
            mutationStrategy = &Algorithms::deMutationStrategy5or10;
            crossoverStrategy = &Algorithms::binomialCrossover;
            break;
    }


    double* popFitness = new double[populationSize];
    double** noisyVectors = new double*[populationSize];
    for(int i=0; i<populationSize; i++){
        popFitness[i] = std::invoke(Algorithms::equation, Algorithms::eq, pop.matrix[i], Algorithms::dimension);
        noisyVectors[i] = new double[Algorithms::dimension];
    }

    for(int i=0; i<generations; i++){

        //create noisy vector
        std::invoke(mutationStrategy, this,  &pop, popFitness, noisyVectors, populationSize, scalingFactor, lambda, &g1);

        //crossover and append population
        std::invoke(crossoverStrategy, this, &pop, popFitness, noisyVectors, populationSize, crossoverFactor, &g1);

    }

    //find best vector in population
    double bestFitness = popFitness[0];
    int bestVector = 0;
    for(int i=1; i< populationSize; i++){
        double newBestFitness = popFitness[i];

        if(newBestFitness < bestFitness){
            bestFitness = newBestFitness;
            bestVector = i;
        }
    }

    //garbage collection
    for(int i=0; i<populationSize; i++){
        delete[] noisyVectors[i];
    }
    delete[] noisyVectors;
    delete[] popFitness;

    for(int i=0; i<Algorithms::dimension; i++){
        bestFit[i] = pop.matrix[bestVector][i];
    }

    return bestFitness;
}


/**
 * performs particle swarm optimization
 * @param particles number of particles
 * @param c1 first learning factor
 * @param c2 second learning factor
 * @param generations
 * @return best fitness found
 */
double Algorithms::particleSwarmOptimization(int particles, double c1, double c2, int generations){

    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::uniform_real_distribution<> velocityDistribution(0,0.5*(Algorithms::upperBound - Algorithms::lowerBound));
    std::uniform_real_distribution<> lD(0,1);
    std::mt19937 g1(seed1);//mersenne_twister_engine

    //allocate memory
    double** particleVelocities = new double*[particles];

    Population pop(particles, Algorithms::dimension, Algorithms::lowerBound, Algorithms::upperBound);
    double* particleFitness = new double[particles];

    double** pBestParticle = new double*[particles];
    double* pBestFitness = new double[particles];

    double* gBestParticle = new double[particles];
    double gBestFitness;

    for(int i=0; i<particles; i++){
        particleVelocities[i] = new double[Algorithms::dimension];
        pBestParticle[i] = new double[Algorithms::dimension];
    }

    //set default values
    for(int i=0; i<particles; i++){
        for(int j=0; j<Algorithms::dimension; j++){
            particleVelocities[i][j] = velocityDistribution(g1);
            pBestParticle[i][j] = pop.matrix[i][j];
        }

        particleFitness[i] = std::invoke(Algorithms::equation, Algorithms::eq, pop.matrix[i], Algorithms::dimension);
        pBestFitness[i] = particleFitness[i];

    }

    //find best fitness
    gBestFitness = pBestFitness[0];
    int bestLocation = 0;
    for(int i=1; i<particles; i++){
        if(pBestFitness[i]< gBestFitness){
            gBestFitness = pBestFitness[i];
            bestLocation = i;
        }
    }

    for(int i=0; i<Algorithms::dimension; i++){
        gBestParticle[i] = pBestParticle[bestLocation][i];
    }

    //particle swarm optimization
    for(int i=0; i<generations; i++){
        for(int j=0; j<particles; j++){
            for(int k=0; k<Algorithms::dimension; k++){
                particleVelocities[j][k] += c1 * lD(g1) * (pBestParticle[j][k] - pop.matrix[j][k]) + c2 * lD(g1) * (gBestParticle[k] - pop.matrix[j][k]);
                pop.matrix[j][k] +=  particleVelocities[j][k];
            }

            forceInBounds(pop.matrix[j]);
            particleFitness[j] = std::invoke(Algorithms::equation, Algorithms::eq, pop.matrix[j], Algorithms::dimension);

            if(particleFitness[j] < pBestFitness[j]){
                pBestFitness[j] = particleFitness[j];

                for(int k=0; k<Algorithms::dimension; k++){
                    pBestParticle[j][k] = pop.matrix[j][k];
                }

                if(particleFitness[j] < gBestFitness){

                    gBestFitness = particleFitness[j];

                    for(int k=0; k<Algorithms::dimension; k++){
                        gBestParticle[k] = pop.matrix[j][k];
                    }

                }

            }

        }
    }

    for(int i=0; i<Algorithms::dimension; i++){
        Algorithms::bestFit[i] = gBestParticle[i];
    }

    //garbage collection
    for(int i=0; i<particles; i++){
        delete[] particleVelocities[i];
        delete[] pBestParticle[i];
    }
    delete[] particleVelocities;
    delete[] pBestParticle;
    delete[] particleFitness;
    delete[] pBestFitness;
    delete[] gBestParticle;

    return gBestFitness;
}













