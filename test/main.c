#include "genalg.h"
#include <stdio.h>
#include <stdint.h>

//linux - unfortunately this doesn't appear to work on my kernel, will need redo-ing
#include <sys/resource.h>
int GetMemUsage(void) {
    struct rusage resourceUsage;
    getrusage(RUSAGE_SELF, &resourceUsage);
    return resourceUsage.ru_idrss;
}

// Code for fitness function copied from usage example!
enum bool checkIfPrime(uint32_t x) {
    if (x <= 3) return true;
    for (uint32_t i = 2; i < x; i++) {
        if (x % i == 0) return false;
    }
    return true;
}

// optimise to find the HCF between two numbers
// this is intentionally poorly optimised to highlight the use-case for multithreaded (performance-intensive fitness function)
double eval(const byte* const genes) {
    // cast and deference 4 bytes per 32b integer
    const uint32_t moduloOperand = 100000;
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
    int result = 0;

    printf("Running unit tests...\n");

    GeneticAlgorithm_t ga = GeneticAlgorithm(
        4000,       // population size
        85,         // % mutation rate
        25,         // % crossover rate
        2,          // ts group size
        8,          // genes
        1200,       // generations
        &eval       // defined fitness evaluation fn
    );

    RunGeneticAlgorithm(&ga, 4, false);
    int memoryUsageFromRun1 = GetMemUsage();
    FreeGeneticAlgorithm(&ga);

    ga = GeneticAlgorithm(
        4000,       // population size
        85,         // % mutation rate
        25,         // % crossover rate
        2,          // ts group size
        8,          // genes
        1200,       // generations
        &eval       // defined fitness evaluation fn
    );

    RunGeneticAlgorithm(&ga, 4, false);
    int memoryUsageFromRun2 = GetMemUsage();

    const uint32_t moduloOperand = 100000;
    // cast and deference 4 bytes per 32b integer
    const uint32_t numbers[2] = { *((uint32_t*)(&(ga.bestSolution.genes[0]))), *((uint32_t*)(&(ga.bestSolution.genes[4])))};
    const double fitness = ga.bestSolution.fitness;

    double res = 0.0;
    if(
        numbers[0] != numbers[1]
        &&
        checkIfPrime(numbers[0] % moduloOperand)
        &&
        checkIfPrime(numbers[1] % moduloOperand)
    ) res = numbers[0] > numbers[1]
        ? numbers[0] - numbers[1]
        : numbers[1] - numbers[0];

    printf("TEST: Fitness value calculated correctly...");
    if (res != fitness) {
        printf("FAIL!\n\tDifference: %f\n\tFitness: %f\n", res, fitness);
        result |= 32;
    } else printf("PASS.\n");

    
    printf("TEST: Ensure no memory leak...");
    // if there is no memory leak, the memory that the OS reserved for Run1 will be reused for Run2, therefore little to no process memory change.
    int changeInMemory = memoryUsageFromRun2 - memoryUsageFromRun1;
    if (changeInMemory > 256) {
        printf("FAIL!\n\tMemory difference: %d\n", changeInMemory);
        result |= 16;
    } else printf("PASS.\n");

    printf("TEST: Optimisation working / fitness value close to 0...");
    if (fitness < -350.0) {
        printf("FAIL!\n\tFitness: %f\n", fitness);
        result |= 64;
    } else printf("PASS.\n");
    
    FreeGeneticAlgorithm(&ga);
    return result;
}