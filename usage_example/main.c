#include "genalg.h"
#include <stdio.h>
#include <stdint.h>

// optimise to minimise the difference between the two numbers defined in the genes
double eval(const byte* const genes) {
    double difference;
    // cast and deference 4 bytes per 32b integer
    const uint32_t numbers[2] = { *((uint32_t*)(&genes[0])), *((uint32_t*)(&genes[4]))}; 

    if (numbers[0] > numbers[1]) {
        difference = (numbers[0] - numbers[1]);
    } else {
        difference = (numbers[1] - numbers[0]);
    }    
    return difference * (-1);
}

int main(void) {
    // init algorithm with:
    GeneticAlgorithm_t ga = GeneticAlgorithm(
        2000,       // 2000 population size
        85000,      // 85,000/100,000 (85%) mutation rate
        25000,      // 25,000/100,000 (25%) crossover rate
        3,          // 3 tournament selection group size
        8,          // 8 genes (1 gene = 1 byte)
        200,        // 200 generations
        &eval       // we point to our defined fitness evaluation "eval"
    );

    RunGeneticAlgorithm(&ga, true); // run with console output=true


    // cast and deference 4 bytes per 32b integer
    const uint32_t numbers[2] = { *((uint32_t*)(&(ga.bestSolution.genes[0]))), *((uint32_t*)(&(ga.bestSolution.genes[4])))}; 
    printf("\nBest solution genes: %u, %u\n\n", numbers[0], numbers[1]);
    FreeGeneticAlgorithm(&ga);

    return 0;
}