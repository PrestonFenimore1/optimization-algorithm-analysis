This program is for comparing the effectiveness of particle swarm optimization and differential evolution for finding the best vector
to minimize a given test function of a provided dimensionality

This program accepts a file in the format

Equation (integer)
Dimensionality (integer)
Population Size (integer)
Search Algorithm (integer)
Lower Bound (double)
Upper Bound (double)
Generations (integer)

as a command line argument

the number to equation list is as follows:

1:Schwefel's function
2:1st De Jong's function
3:Rosenbrock
4:Rastrigin
5:Griewangk
6:Sine Envelope Sine Wave
7:Stretched V Sine Wave
8:Ackley's One
9:Ackley's Two
10:Egg Holder

the number to search algorithm is as follows

1:Differential Evolution
2:Particle Swarm Optimization

additional parmeters are also needed depending on the search algorithm

for Differential Evolution

Crossover Factor (double)
Scaling Factor (double)
Second Scaling Factor (double)(set as 0 when unneeded)
Evolution Strategy (integer) (see details below)

Accepted strategies are as follows

1:best/1/exp
2:rand/1/exp
3:rand-to-best/1/exp
4:best/2/exp
5:rand/2/exp
6:best/1/bin
7:rand/1/bin
8:rand-to-best/1/bin
9:best/2/bin
10:rand/2/bin

where bin is binomial crossover, exp is exponential crossover, and x and y in the prior data in the form x/y/z is 
the vector to be perturbed and the number of difference vectors considered for perturbation

for Particle Swarm Optimization

First Learning Factor (double)
Second Learning Factor (double)


The output is a .csv file of the best vector found and its fitness value.
Example input file "sampleInput1.txt" and "sampleInput2.txt" have been provided in the source folder.
They are an example for Differential Evolution and Particle Swarm Optimization respectively

