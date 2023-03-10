#include "genalg.h"
#include <stdio.h>
#include <stdint.h>

#ifdef _WIN32
    #include <windows.h>
    #include <psapi.h>

    int GetMemUsage(void) {
        HANDLE currentProcPseudohandle = GetCurrentProcess();
        PROCESS_MEMORY_COUNTERS_EX procMemoryCtrs;

        if (GetProcessMemoryInfo( currentProcPseudohandle, &procMemoryCtrs, sizeof(procMemoryCtrs)))
        {
            return procMemoryCtrs.PrivateUsage;
        }

        return 0;
    }
#else
	//linux - unfortunately this doesn't appear to work on my kernel, so I'll abandon this for now.
    #include <sys/resource.h>

    int GetMemUsage(void) {
        struct rusage resourceUsage;
        getrusage(RUSAGE_SELF, &resourceUsage);
        return resourceUsage.ru_idrss;
    }
#endif


// much code copied from usage_example
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
    int result = 0;

    GeneticAlgorithm_t ga = GeneticAlgorithm(
        8000,       // population size
        85,         // % mutation rate
        25,         // % crossover rate
        2,          // ts group size
        8,          // genes
        2000,        // generations
        &eval       // defined fitness evaluation fn
    );

    RunGeneticAlgorithm(&ga, false);
    int memoryUsageFromRun1 = GetMemUsage();
    FreeGeneticAlgorithm(&ga);

    ga = GeneticAlgorithm(
        8000,       // population size
        85,         // % mutation rate
        25,         // % crossover rate
        2,          // ts group size
        8,          // genes
        2000,        // generations
        &eval       // defined fitness evaluation fn
    );

    RunGeneticAlgorithm(&ga, false);
    int memoryUsageFromRun2 = GetMemUsage();

    // cast and deference 4 bytes per 32b integer
    const uint32_t numbers[2] = { *((uint32_t*)(&(ga.bestSolution.genes[0]))), *((uint32_t*)(&(ga.bestSolution.genes[4])))};
    const double fitness = ga.bestSolution.fitness;
    double difference;
    
    // TEST: Fitness value matches return
    if (numbers[0] > numbers[1]) {
         difference = (numbers[0] - numbers[1]);
    } else {
         difference = (numbers[1] - numbers[0]);
    }
    difference *= (-1.0);
    if (difference != fitness) {
        printf("Fitness == fitnessEval test failed! Difference: %f, Fitness: %f\n", difference, fitness);
        result |= 32;
    }

    // TEST: Memory leak
    // if there is no memory leak, the memory that the OS reserved for Run1 will be reused for Run2, therefore little to no process memory change.
    int changeInMemory = memoryUsageFromRun2 - memoryUsageFromRun1;
    if (changeInMemory > 256) {
        printf("Memory leak test failed! Memory difference: %d\n", changeInMemory);
        result |= 16;
    }

    // TEST: Fitness is close to 0 (i.e. shows optimisation has occurred)
    if (fitness < -350.0) {
        printf("Fitness threshold test failed! Fitness: %f\n", fitness);
        result |= 64;
    }
    
    FreeGeneticAlgorithm(&ga);
    return result;
}