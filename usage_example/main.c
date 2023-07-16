#include "genalg.h"
#include <stdio.h>
#include <stdint.h>

enum bool checkIfPrime(uint32_t x) {
    if (x <= 3) return true;
    for (uint32_t i = 2; i < x; i++) {
        if (x % i == 0) return false;
    }
    return true;
}

// optimise to find the HCF between two numbers and difference
// this is intentionally poorly optimised to highlight the use-case for multithreaded (performance-intensive fitness function)
// runs approx. 3x faster on 4 threads than single-threaded, HOWEVER if using a very simple, low CPU-time fitness function, the difference will be less
// with some fitness functions, running single-threaded will actually be faster
double eval(const byte* const genes) {
    // cast and deference 4 bytes per 32b integer
    const uint32_t moduloOperand = 50000;
    const uint32_t numbers[2] = { *((uint32_t*)(genes)), *((uint32_t*)(genes + 4)) };

    if(
        numbers[0] != numbers[1]
        &&
        checkIfPrime(numbers[0] % moduloOperand)
        &&
        checkIfPrime(numbers[1] % moduloOperand)
    ) return numbers[0] > numbers[1]
        ? numbers[0] - numbers[1]
        : numbers[1] - numbers[0];
    else return 0.0;
}

int main(void) {
    // init algorithm with:
    GeneticAlgorithm_t ga = GeneticAlgorithm(
        20000,      // 20000 population size
        40,         // 40/100 (40%) mutation rate
        30,         // 30/100 (30%) crossover rate
        2,          // 2 tournament selection group size
        8,          // 8 genes (1 gene = 1 byte)
        200,        // 200 generations
        &eval       // we point to our defined fitness evaluation "eval"
    );

    // run with 4 threads, and console output=true
    // set the second argument to 1, to run SINGLE-THREADED mode
    RunGeneticAlgorithm(&ga, 4, true);



    // cast and deference 4 bytes per 32b integer
    const uint32_t numbers[2] = { *((uint32_t*)(&(ga.bestSolution.genes[0]))), *((uint32_t*)(&(ga.bestSolution.genes[4])))}; 
    printf("\nBest solution fitness / genes:\n\t%f / %u, %u\n\n", ga.bestSolution.fitness, numbers[0], numbers[1]);
    FreeGeneticAlgorithm(&ga);

    return 0;
}