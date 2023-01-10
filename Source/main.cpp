#include <iostream>
#include <fstream>
#include "Algorithms.h"

/**
 * saves equation name to equationName based on equation number
 * @param equationName pointer to equationName
 * @param equationNumber number representing equation
 */
void getEquationName(std::string* equationName, int equationNumber){

    switch(equationNumber){
        case 1:
            *equationName = "Schwefel's";
            break;
        case 2:
            *equationName = "1st De Jong's";
            break;
        case 3:
            *equationName = "Rosenbrock";
            break;
        case 4:
            *equationName = "Rastrigin";
            break;
        case 5:
            *equationName = "Griewangk";
            break;
        case 6:
            *equationName = "Sine Envelope Sine Wave";
            break;
        case 7:
            *equationName = "Stretched V Sine Wave";
            break;
        case 8:
            *equationName = "Ackley's One";
            break;
        case 9:
            *equationName = "Ackley's Two";
            break;
        case 10:
            *equationName = "Egg Holder";
            break;
    }

}

/**
 * gets strategy name from number
 * @param strategyToPrint
 * @param strategyNum
 */
void getStrategy(std::string* strategyToPrint, int strategyNum){

    switch(strategyNum){
        case 1:
            *strategyToPrint = "DE/best/1/exp";
            break;
        case 2:
            *strategyToPrint = "DE/rand/1/exp";
            break;
        case 3:
            *strategyToPrint = "DE/rand-to-best/1/exp";
            break;
        case 4:
            *strategyToPrint = "DE/best/2/exp";
            break;
        case 5:
            *strategyToPrint = "DE/rand/2/exp";
            break;
        case 6:
            *strategyToPrint = "DE/best/1/bin";
            break;
        case 7:
            *strategyToPrint = "DE/rand/1/bin";
            break;
        case 8:
            *strategyToPrint = "DE/rand-to-best/1/bin";
            break;
        case 9:
            *strategyToPrint = "DE/best/2/bin";
            break;
        case 10:
            *strategyToPrint = "DE/rand/2/bin";
            break;

    }

}

/**
 * performs differential evolution and outputs result
 * @param equationNumber
 * @param dimension
 * @param lowerBound
 * @param upperBound
 * @param populationSize
 * @param crossoverFactor
 * @param scalingFactor
 * @param lambda
 * @param generations
 * @param strategy
 */
void getDEData(int equationNumber, int dimension, double lowerBound, double upperBound, int populationSize, double crossoverFactor,
    double scalingFactor, double lambda, int generations, int strategy){

    Algorithms alg(equationNumber, dimension, lowerBound, upperBound);
    double fitness = alg.differentialEvolution(populationSize, crossoverFactor, scalingFactor, lambda, generations, strategy);

    std::string dimensionToPrint = std::to_string(dimension);
    std::string equationName;
    getEquationName(&equationName, equationNumber);
    std::string strategyToPrint;
    getStrategy(&strategyToPrint, strategy);

    std::string vector = std::to_string(alg.bestFit[0]);

    for(int i=1; i<dimension; i++){
        vector.append("," + std::to_string(alg.bestFit[i]));
    }

    std::string output = dimensionToPrint + "," + equationName + "," +
        strategyToPrint + "," + std::to_string(fitness) + "," + vector;


    std::ofstream outputFile("output.csv", std::ios_base::app);

    outputFile << output + "\n";

    outputFile.close();

}

/**
 * performs particle swarm optimization and outputs result
 * @param equationNumber
 * @param dimension
 * @param lowerBound
 * @param upperBound
 * @param populationSize
 * @param firstLearningFactor
 * @param secondLearningFactor
 * @param generations
 */
void getPSData(int equationNumber, int dimension, double lowerBound, double upperBound,
    int populationSize, double firstLearningFactor, double secondLearningFactor, int generations){

    Algorithms alg(equationNumber, dimension, lowerBound, upperBound);
    double fitness = alg.particleSwarmOptimization(populationSize, firstLearningFactor, secondLearningFactor, generations);

    std::string dimensionToPrint = std::to_string(dimension);
    std::string equationName;
    getEquationName(&equationName, equationNumber);

    std::string vector = std::to_string(alg.bestFit[0]);

    for(int i=1; i<dimension; i++){
        vector.append("," + std::to_string(alg.bestFit[i]));
    }

    std::string output = dimensionToPrint + "," + equationName + "," +
        "Particle Swarm" + "," + std::to_string(fitness) + "," + vector;

    std::ofstream outputFile("output.csv", std::ios_base::app);

    outputFile << output + "\n";

    outputFile.close();

}

/**
 * accepts input file and performs provided algorithm with given parameters
 * then saves to file "output.csv"
 * @param argc
 * @param argv filename
 * @return exit code
 */
int run(int argc, char* argv[]){

    if(argc != 2){
        std::cout << "incorrect number of operands" << std::endl;
        return 1;
    }

    std::ifstream inputFile(argv[1]);

    char equationInput[256];
    char dimensionalityInput[256];
    char populationSizeInput[256];
    char algorithmInput[256];
    char lowerBoundInput[256];
    char upperBoundInput[256];
    char timesToRunInput[256];

    inputFile.getline(equationInput,256);
    inputFile.getline(dimensionalityInput,256);
    inputFile.getline(populationSizeInput,256);
    inputFile.getline(algorithmInput,256);
    inputFile.getline(lowerBoundInput,256);
    inputFile.getline(upperBoundInput,256);
    inputFile.getline(timesToRunInput,256);

    int equation = std::stoi(equationInput);
    int dimensionality = std::stoi(dimensionalityInput);
    int populationSize = std::stoi(populationSizeInput);
    int algorithm = std::stoi(algorithmInput);
    double lowerBound = std::stod(lowerBoundInput);
    double upperBound = std::stod(upperBoundInput);
    int timesToRun = std::stoi(timesToRunInput);

    switch(algorithm){
        case 1:
            {
                char crossoverFactorInput[256];
                char scalingFactorInput[256];
                char secondScalingFactorInput[256];
                char strategyInput[256];

                inputFile.getline(crossoverFactorInput, 256);
                inputFile.getline(scalingFactorInput, 256);
                inputFile.getline(secondScalingFactorInput, 256);
                inputFile.getline(strategyInput, 256);

                double crossoverFactor = std::stod(crossoverFactorInput);
                double scalingFactor = std::stod(scalingFactorInput);
                double secondScalingFactor = std::stod(secondScalingFactorInput);
                int strategy = std::stoi(strategyInput);
                inputFile.close();

                getDEData(equation, dimensionality, lowerBound, upperBound, populationSize, crossoverFactor,
                    scalingFactor, secondScalingFactor, timesToRun, strategy);
            }
            break;
        case 2:
            {
                char firstLearningFactorInput[256];
                char secondLearningFactorInput[256];

                inputFile.getline(firstLearningFactorInput, 256);
                inputFile.getline(secondLearningFactorInput, 256);

                double firstLearningFactor = std::stod(firstLearningFactorInput);
                double secondLearningFactor = std::stod(secondLearningFactorInput);
                inputFile.close();

                getPSData(equation, dimensionality, lowerBound, upperBound, populationSize, firstLearningFactor,
                          secondLearningFactor, timesToRun);
            }
            break;
        default:
            std::cout << "algorithm must be 1 or 2" << std::endl;
            inputFile.close();
            return 1;
    }

    return 0;

}

/**
 *
 * @param argc
 * @param argv filename
 * @return exit code
 */
int main(int argc, char* argv[]) {

    return run( argc, argv);
}
