#include "genalg.h"
#include <stdio.h>
#include <stdint.h>

// optimise the minimise the difference between the two numbers defined in each gene
double eval(const uint32_t* const genes) {
    double difference;
    if (genes[0] > genes[1]) {
        difference = (genes[0] - genes[1]);
    } else {
        difference = (genes[1] - genes[0]);
    }    
    return difference * (-1);
}

int main() {
    // init algorithm with:
    GeneticAlgorithm_t ga = GeneticAlgorithm(
        200,        // 200 population size
        70000,      // 70,000/100,000 (70%) mutation rate
        35000,      // 35,000/100,000 (35%) crossover rate
        4,          // 4 tournament selection group size
        2,          // 2 genes
        240,        // 240 generations
        &eval       // we point to our defined fitness evaluation "eval"
    );

    RunGeneticAlgorithm(&ga, true); // run with console output=true
    printf("\nBest solution genes: %u, %u\n\n", ga.bestSolution.genes[0], ga.bestSolution.genes[1]);
    FreeGeneticAlgorithm(&ga);

    return 0;
}