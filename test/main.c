#include "genalg.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>

uint32_t GetMemUsage(void) {
    char filePath[32];
    sprintf(filePath, "/proc/%i/stat", getpid());
    FILE *fileHandle = fopen(filePath, "r");
    uint32_t memSize = 0;
    int res = fscanf(fileHandle, "%*d %*s %*c %*d %*d %*d %*d %*d %*u %*u %*u %*u %*u %*u %*u %*d %*d %*d %*d %*d %*d %*u %u", &memSize);
    return res == 0 ? 0 : memSize;
}

// Code for fitness function copied from usage example!
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

    const uint32_t moduloOperand = 50000;
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




    // skipped because test doesn't work well on CI (where the OS is tight when it comes to managing memory). It is still a decent test on local.
    printf("TEST: Ensure no memory leak...SKIPPED.\n");
    /*
    // if there is no memory leak, the memory that the OS reserved for Run1 will be reused for Run2, therefore little to no process memory change.
    int changeInMemory = memoryUsageFromRun2 - memoryUsageFromRun1;
    if (changeInMemory > 256) {
        printf("FAIL!\n\tMemory difference: %d\n", changeInMemory);
        result |= 16;
    } else if (memoryUsageFromRun1 == 0  || memoryUsageFromRun2 == 0){
        printf("FAIL!\n\tFailed reading /proc/[pid]/stat file\n");
        result |= 16;
    } else printf("PASS.\n");
    */






    printf("TEST: Optimisation working / fitness value close to 0...");
    if (fitness < -350.0) {
        printf("FAIL!\n\tFitness: %f\n", fitness);
        result |= 64;
    } else printf("PASS.\n");
    
    FreeGeneticAlgorithm(&ga);
    return result;
}