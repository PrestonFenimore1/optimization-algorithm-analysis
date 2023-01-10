#include "Population.h"

/**
 * creates a Population object with a generated matrix as a class variable
 * @param experiments number of rows in generated population matrix
 * @param dimensionality  number of columns in generated population matrix
 * @param lowerBound lower bound for values in population matrix
 * @param upperBound upper bound bound for values in population matrix
 */
Population::Population(int experiments, int dimensionality, double lowerBound, double upperBound){

    Population::rows = experiments;
    Population::columns = dimensionality;

    Population::matrix = new double*[rows];
    Population::fitness = new double[rows];

    for(int i=0; i< rows; i++){
        matrix[i] = new double[columns];
    }

    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 g1 (seed1);//mersenne_twister_engine

    std::uniform_real_distribution<> d(lowerBound,upperBound);

    for(int i=0; i< rows; i++){
        for(int j=0; j< columns; j++){

            matrix[i][j] = d(g1);

        }
    }

}

Population::~Population(){
    for(int i=0; i< rows; i++){
        delete[] matrix[i];
    }
    delete[] matrix;
    delete[] fitness;
}

/**
 * runs vectors in population matrix through given equation and stores fitness vector in Population::fitness
 * @param type which equation will be used for experiment
 */
void Population::runExperiment(int type){

    Equations e;

    switch (type){
        case 1:
            for(int i=0; i<rows; i++){
                Population::fitness[i] = e.schwefelsfunction(matrix[i],columns);
            }
            break;
        case 2:
            for(int i=0; i<rows; i++){
                Population::fitness[i] = e.firstDeJongsfunction(matrix[i],columns);
            }
            break;
        case 3:
            for(int i=0; i<rows; i++){
                Population::fitness[i] = e.rosenbrock(matrix[i],columns);
            }
            break;
        case 4:
            for(int i=0; i<rows; i++){
                Population::fitness[i] = e.rastrigin(matrix[i],columns);
            }
            break;
        case 5:
            for(int i=0; i<rows; i++){
                Population::fitness[i] = e.griewangk(matrix[i],columns);
            }
            break;
        case 6:
            for(int i=0; i<rows; i++){
                Population::fitness[i] = e.sineEnvelopeSineWave(matrix[i],columns);
            }
            break;
        case 7:
            for(int i=0; i<rows; i++){
                Population::fitness[i] = e.stretchedVSineWave(matrix[i],columns);
            }
            break;
        case 8:
            for(int i=0; i<rows; i++){
                Population::fitness[i] = e.ackleysOne(matrix[i],columns);
            }
            break;
        case 9:
            for(int i=0; i<rows; i++){
                Population::fitness[i] = e.ackleysTwo(matrix[i],columns);
            }
            break;
        case 10:
            for(int i=0; i<rows; i++){
                Population::fitness[i] = e.eggHolder(matrix[i],columns);
            }
            break;
        default:
            throw std::invalid_argument( "equationNumber must be 1-10" );
    }

}

