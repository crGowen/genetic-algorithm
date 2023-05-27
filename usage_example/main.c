#include "genalg.h"
#include <stdio.h>
#include <stdint.h>

// optimise to minimise the difference between the two numbers defined in the genes
double eval(const byte* const genes) {
    double difference;
    // cast and deference 4 bytes per 32b integer
    const uint32_t numbers[2] = { *((uint32_t*)(genes)), *((uint32_t*)(genes + 4)) };

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
        20000,      // 20000 population size
        60,         // 60/100 (60%) mutation rate
        25,         // 25/100 (25%) crossover rate
        2,          // 2 tournament selection group size
        8,          // 8 genes (1 gene = 1 byte)
        200,        // 200 generations
        &eval       // we point to our defined fitness evaluation "eval"
    );

    RunGeneticAlgorithm(&ga, 4, true); // run with 4 threads, and console output=true


    // cast and deference 4 bytes per 32b integer
    const uint32_t numbers[2] = { *((uint32_t*)(&(ga.bestSolution.genes[0]))), *((uint32_t*)(&(ga.bestSolution.genes[4])))}; 
    printf("\nBest solution genes: %u, %u\n\n", numbers[0], numbers[1]);
    FreeGeneticAlgorithm(&ga);

    return 0;
}