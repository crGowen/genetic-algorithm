#include "genalg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

#define CACHE_LINE_SIZE 64

typedef struct {
    double (*fitnessEvaluationFn)(const byte * const genes);
    byte* populationGenes;
    byte* nextGenGenes;
    double* populationFitness;
    pthread_barrier_t* barrier;
    uint32_t totalPopSize;
    uint16_t numberOfGenerations;
    uint16_t groupSize;
    uint16_t numberOfGenes;
    uint8_t mutationRate;
    uint8_t crossoverRate;
    enum bool enablePrintOutput;
} CommonGaThreadArgs_t;

typedef struct {
    CommonGaThreadArgs_t* commonArgs;
    uint32_t start;
    uint32_t popSize;
    uint8_t threadId;
} GeneticAlgorithmThread_t;

uint32_t GenerateRandomU32(void);
byte GenerateRandomByte(void);
void* RunGenAlgThread(void* rawArgs);
void FreeGeneticAlgorithm(GeneticAlgorithm_t* const ga);

// end of function declarations




void RunGeneticAlgorithm(GeneticAlgorithm_t* const ga, const uint8_t numThreads, const enum bool enablePrintOutput) {
    if (ga->block) {
		printf("\nGenetic algorithm not initialised properly. Will not run.");
		return;
	}

    const uint32_t popThreadFactor = CACHE_LINE_SIZE * numThreads;

    ga->popSize = ga->popSize + popThreadFactor - (ga->popSize % popThreadFactor);
    const uint64_t popGenesMemSize = ga->popSize * ga->numberOfGenes * sizeof(byte);

    byte * populationGenes = (byte*) malloc(popGenesMemSize);
    byte * nextGenGenes = (byte*) malloc(popGenesMemSize);

    ga->populationFitness = (double*) malloc(ga->popSize * sizeof(double));


    for (uint64_t i = 0; i < popGenesMemSize; i++) {
        populationGenes[i] = GenerateRandomByte();
        nextGenGenes[i] = GenerateRandomByte();
    }

    ga->bestSolution.fitness = 0.0;

    GeneticAlgorithmThread_t* threadArgs = (GeneticAlgorithmThread_t*) malloc(numThreads * sizeof(GeneticAlgorithmThread_t));

    const uint32_t popPerThread = ga->popSize / numThreads;
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, numThreads);

    CommonGaThreadArgs_t commonArgs = (CommonGaThreadArgs_t){
        .totalPopSize          = ga->popSize,
        .fitnessEvaluationFn   = ga->fitnessEvaluationFn,
        .barrier               = &barrier,
        .groupSize             = ga->groupSize,
        .numberOfGenes         = ga->numberOfGenes,
        .numberOfGenerations   = ga->numberOfGenerations,
        .mutationRate          = ga->mutationRate,
        .crossoverRate         = ga->crossoverRate,
        .populationFitness     = ga->populationFitness,
        .enablePrintOutput     = enablePrintOutput,

        .populationGenes = populationGenes,
        .nextGenGenes = nextGenGenes
    };

    for (int i = 0; i < numThreads; i++) {
        threadArgs[i].commonArgs            = &commonArgs;
        threadArgs[i].threadId              = i;
        threadArgs[i].start                 = i * popPerThread;
        threadArgs[i].popSize               = popPerThread;
    }

    pthread_t* threads = (pthread_t*) malloc((numThreads - 1) * sizeof(pthread_t));

    for (int i = 0; i < numThreads - 1; i++) {
        pthread_create(threads + i, NULL, RunGenAlgThread, (void*)(threadArgs + i + 1));
    }

    RunGenAlgThread(threadArgs);

    for (int i = 0; i < numThreads - 1; i++) {
        pthread_join(threads[i], NULL);
    }

    pthread_barrier_destroy(&barrier);

    ga->populationGenes = commonArgs.populationGenes;

    // final fitness evaluation
    enum bool bestFitnessIsAssigned = false;
    for (uint32_t i = 0; i < ga->popSize; i++) {
        if (!bestFitnessIsAssigned || ga->populationFitness[i] > ga->bestSolution.fitness) {
            bestFitnessIsAssigned = true;
            ga->bestSolution.fitness = ga->populationFitness[i];
            ga->bestSolution.genes = ga->populationGenes + (i * ga->numberOfGenes);
        }
    }

    free(commonArgs.nextGenGenes);
    free(threadArgs);
    free(threads);
}

void* RunGenAlgThread(void* rawArgs) {
    GeneticAlgorithmThread_t* gaThread = (GeneticAlgorithmThread_t*)rawArgs;

    const uint16_t gens = gaThread->commonArgs->numberOfGenerations;
    const uint32_t byteOffset = (gaThread->start * gaThread->commonArgs->numberOfGenes);
    const uint32_t fitnessOffset = gaThread->start;
    const uint8_t threadId = gaThread->threadId;
    const uint32_t popSize = gaThread->popSize;
    const uint32_t totalPopSize = gaThread->commonArgs->totalPopSize;
    const uint16_t numberOfGenes = gaThread->commonArgs->numberOfGenes;
    const uint16_t groupSize = gaThread->commonArgs->groupSize;
    const uint8_t crossoverRate = gaThread->commonArgs->crossoverRate;
    const uint8_t mutationRate = gaThread->commonArgs->mutationRate;
    const uint8_t bitsPerGene = 8 * sizeof(byte);
    const enum bool enablePrintOutput = gaThread->commonArgs->enablePrintOutput;

    uint32_t* contenders = (uint32_t*) malloc(groupSize * sizeof(uint32_t));

    uint32_t bitPosition, crossMask, crossIndex, mutationMask, mutationIndex;

    const double (*fitnessEvaluationFn)(const byte * const genes) = gaThread->commonArgs->fitnessEvaluationFn;
    double const * populationFitness = gaThread->commonArgs->populationFitness;

    for (uint32_t currentGen = 0; currentGen < gens; currentGen++) {
        byte*  populationGenes = gaThread->commonArgs->populationGenes;
        byte*  threadPop = populationGenes + byteOffset;
        byte*  threadNext = gaThread->commonArgs->nextGenGenes + byteOffset;
        double* threadFitness = populationFitness + fitnessOffset;
        double bestFitnessForThread = 0.0;
        enum bool bestFitnessIsAssigned = false;

        // fitness evaluation
        for (uint32_t i = 0; i < popSize; i++) {
            threadFitness[i] = fitnessEvaluationFn(
                threadPop + (i * numberOfGenes)
            );
            if (!bestFitnessIsAssigned || threadFitness[i] > bestFitnessForThread) {
                bestFitnessIsAssigned = true;
                bestFitnessForThread = threadFitness[i];
            }
        }

        pthread_barrier_wait(gaThread->commonArgs->barrier);

        // tournament selection - needs DRYing
        const byte *mother, *father;
        uint32_t motherIndex;
        for (uint32_t childIndex = 0; childIndex < popSize; childIndex++ ) {
            // MOTHER
            // reset contenders
            for(uint32_t i = 0; i < groupSize; i++) {
                contenders[i] = 0;
            }

            // find mother contenders
            for (uint32_t i = 0; i < groupSize; i++) {
                enum bool isUniqueContender = false;
                uint32_t randomSolution;
                while (isUniqueContender == false) {
                    isUniqueContender = true;
                    randomSolution = GenerateRandomU32() % totalPopSize;
                    for (uint32_t j = 0; isUniqueContender == true && j < groupSize; j++) {
                        if (contenders[j] == randomSolution) isUniqueContender = false;
                    }
                }
                contenders[i] = randomSolution;
            }
            // compare contenders and find the best
            uint32_t bestContender = contenders[0];
            for (uint32_t i = 0; i < groupSize; i++) {
                if (populationFitness[contenders[i]] > populationFitness[bestContender]) bestContender = contenders[i];
            }
            motherIndex = bestContender;
            mother = populationGenes + (bestContender * numberOfGenes);

            // FATHER
            // reset contenders
            for(uint32_t i = 0; i < groupSize; i++) {
                contenders[i] = 0;
            }

            // find father contenders
            for (uint32_t i = 0; i < groupSize; i++) {
                enum bool isUniqueContender = false;
                uint32_t randomSolution;
                while (isUniqueContender == false) {
                    randomSolution = GenerateRandomU32() % totalPopSize;
                    isUniqueContender = (randomSolution != motherIndex);
                    for (uint32_t j = 0; isUniqueContender == true && j < groupSize; j++) {
                        if (contenders[j] == randomSolution) isUniqueContender = false;
                    }
                }
                contenders[i] = randomSolution;
            }
            // compare contenders and find the best
            bestContender = contenders[0];
            for (uint32_t i = 0; i < groupSize; i++) {
                if (populationFitness[contenders[i]] > populationFitness[bestContender]) bestContender = contenders[i];
            }
            father = populationGenes + (bestContender * numberOfGenes);

            // prepare to generate child
            crossIndex = numberOfGenes + 1;
            mutationIndex = numberOfGenes + 1;

            // prepare crossover
            if (GenerateRandomU32() % 100 < crossoverRate) {
                // N between 1 and bitsPerGene*numberOfGenes - 1
                bitPosition = 1 + (GenerateRandomU32() % (bitsPerGene * numberOfGenes - 1));
                // index of the byte array in which the bit must be
                crossIndex = bitPosition / bitsPerGene;
                // get bitmask which is all 0s on one side and all 1s on the other side of the bit position
                crossMask = (1 << (bitPosition % bitsPerGene)) - 1;
            }

            // prepare mutation
            if (GenerateRandomU32() % 100 < mutationRate) {
                // N between 1 and bitsPerGene*numberOfGenes - 1
                bitPosition = 1 + (GenerateRandomU32() % (bitsPerGene * numberOfGenes - 1));
                // index of the byte array in which the bit must be
                mutationIndex = bitPosition / bitsPerGene;
                // get bitmask which is all 0s on one side and all 1s on the other side of the bit position
                mutationMask = (1 << (bitPosition % bitsPerGene)) - 1;
            }

            // Generate new child
            byte* child = threadNext + (childIndex * numberOfGenes);
            for (uint32_t i = 0; i < numberOfGenes; i++) {
                // Crossover (where mother and father's DNA meets)
                if (i < crossIndex) child[i] = mother[i];
                else if (i > crossIndex) child[i] = father[i];
                else child[i] = (mother[i] & crossMask) + (father[i] & ~crossMask);

                // Mutation (where a single bit is flipped)
                if (i == mutationIndex) child[i] = (child[i] & ~mutationMask) + ((~child[i]) & mutationMask);
            }
        }

        pthread_barrier_wait(gaThread->commonArgs->barrier);

        if (threadId == 0) {
            //control thread
            if (enablePrintOutput) printf("Generation %u completed. Fitness indication = %f\n", currentGen + 1, bestFitnessForThread);
            byte* temp = gaThread->commonArgs->populationGenes;
            gaThread->commonArgs->populationGenes = gaThread->commonArgs->nextGenGenes;
            gaThread->commonArgs->nextGenGenes = temp;
        }

        pthread_barrier_wait(gaThread->commonArgs->barrier);
    }

    if (threadId == 0) {
        if (enablePrintOutput) printf("All generations completed\n");
        byte* temp = gaThread->commonArgs->populationGenes;
        gaThread->commonArgs->populationGenes = gaThread->commonArgs->nextGenGenes;
        gaThread->commonArgs->nextGenGenes = temp;
    }

    pthread_barrier_wait(gaThread->commonArgs->barrier);

    free(contenders);

    return NULL;
}

GeneticAlgorithm_t GeneticAlgorithm(
    const uint32_t popSize,
    const uint8_t mutationRate,
    const uint8_t crossoverRate,
    const uint16_t groupSize,
    const uint16_t numberOfGenes,
    const uint16_t numberOfGenerations,
    double (* const fitnessEvaluationFn)(const byte* const genes)
) {
    GeneticAlgorithm_t ga = (GeneticAlgorithm_t){
        .fitnessEvaluationFn = fitnessEvaluationFn,
        .popSize = popSize,
        .mutationRate = mutationRate,
        .crossoverRate = crossoverRate,
        .groupSize = groupSize,
        .numberOfGenes = numberOfGenes,
        .numberOfGenerations = numberOfGenerations,
        .block = false,
        .populationFitness = NULL,
        .populationGenes = NULL,
        .bestSolution = (BestSolution_t){
            .fitness = 0.0,
            .genes = NULL
        }
    };

    if (groupSize > popSize) {
        printf("\nGroup size cannot be greater population size!");
        ga.block = true;
    }

    srand((unsigned)time(NULL));
    return ga;
}

void FreeGeneticAlgorithm(GeneticAlgorithm_t* const ga) {
    free(ga->populationGenes);
    free(ga->populationFitness);
}

uint32_t GenerateRandomU32(void) {
    // existing randomisation functions weren't being random enough, not even <random>
	// so this function should be a little more random

	uint32_t u32RandomNumber = 0;

	for (uint32_t i = 0; i < 2; i++) {
		const uint32_t u16IntermediateRandom = rand() % 65536;
		const uint32_t u8Offset = i * 16;
		u32RandomNumber += (u16IntermediateRandom << u8Offset);
	}

	return u32RandomNumber;
}

byte GenerateRandomByte(void) {
	return rand() % 256;
}