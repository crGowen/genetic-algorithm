#include "genalg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// GaSolution constructors/destructors
void FreeGaSolution(GaSolution_t* const solution) {
    if (solution->genes != NULL) {
        free(solution->genes);
    }
}

GaSolution_t CopyGaSolution(const GaSolution_t* const from) {
    GaSolution_t to = *from;

    to.genes = (byte*) malloc(to.numberOfGenes * sizeof(byte));
    memcpy(to.genes, from->genes, to.numberOfGenes * sizeof(byte));
    
    return to;
}

GaSolution_t GaSolution(const uint16_t numberOfGenes, const enum bool isRandomGenes) {
    GaSolution_t solution = (GaSolution_t){.numberOfGenes = numberOfGenes, .fitness = 0.0};

    solution.genes = (byte*) malloc(numberOfGenes * sizeof(byte));

    if (isRandomGenes == true) {
        for (uint64_t i = 0; i < numberOfGenes; i++) {
            solution.genes[i] = GenerateRandomByte();
        }
    }

    return solution;
}

// GeneticAlgorithm constructors/destructors
void FreePopulation(GeneticAlgorithm_t* const ga) {
    if (ga->population != NULL) {
        for (uint64_t i = 0; i < ga->popSize; i++) {
            FreeGaSolution(&(ga->population[i]));
        }
        free(ga->population);
    }
}

void FreeGeneticAlgorithm(GeneticAlgorithm_t* const ga) {
    FreePopulation(ga);
    FreeGaSolution(&(ga->bestSolution));
}

GeneticAlgorithm_t GeneticAlgorithm(
    const uint32_t popSize,
    const uint32_t mutationRate,
    const uint32_t crossoverRate,
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
        .block = false
    };

    if (groupSize > popSize) {
        printf("\nGroup size cannot be greater population size!");
        ga.block = true;
    }

    srand((unsigned)time(NULL));

    return ga;
}



// All other functions
GaSolution_t GenerateChild(const GeneticAlgorithm_t* const ga, const GaSolution_t* const mother, const GaSolution_t* const father) {
    uint32_t bitPosition;
    uint32_t crossMask = 0;
    uint32_t crossIndex = ga->numberOfGenes + 1;
    uint32_t mutateMask = 0;
    uint32_t mutateIndex = ga->numberOfGenes + 1;
    const byte bitsPerGene = 8;

    if (GenerateRandomU32() % 100000 < ga->crossoverRate) {
        // crossover
        bitPosition = 1 + (GenerateRandomU32() % (bitsPerGene * ga->numberOfGenes - 1)); // N between 1 and bitsPerGene*NumberOfGenes - 1 (i.e. bit position between first and last bit of u32 array)
        crossIndex = bitPosition / bitsPerGene; // get the index of the byte array in which the bit must be
		crossMask = (1 << (bitPosition % bitsPerGene)) - 1; // get bitmask which will be all zeroes on 1 side and all 1s of the otherside, of the bit position
    }

    if (GenerateRandomU32() % 100000 < ga->mutationRate) {
        //mutation
		bitPosition = 1 + (GenerateRandomU32() % (bitsPerGene * ga->numberOfGenes - 1)); // N between 1 and bitsPerGene*NumberOfGenes - 1 (i.e. bit position between first and last bit of u32 array)

		mutateIndex = bitPosition / bitsPerGene; // get the index of the byte array in which the bit must be
		mutateMask = (1 << (bitPosition % (bitsPerGene - 1))); // get bitmask which will be all zeroes on 1 side and all 1s of the otherside, of the bit position
    }

    GaSolution_t solution = GaSolution(ga->numberOfGenes, false);

    for (uint32_t i = 0; i < ga->numberOfGenes; i++) {
		// crossover (aka where the mother's DNA meets the father's ... no, not like that)
		if (i < crossIndex)
			solution.genes[i] = mother->genes[i];
		else if (i > crossIndex)
			solution.genes[i] = father->genes[i];
		else
			solution.genes[i] = (mother->genes[i] & crossMask) + (father->genes[i] & ~crossMask);

		// mutation (aka flip a single bit)
		if (i == mutateIndex)
			solution.genes[i] = (solution.genes[i] & ~mutateMask) + (~(solution.genes[i]) & mutateMask);
	}

	return solution;
}

void CreatePopulation(GeneticAlgorithm_t* const ga) {
    ga->population = (GaSolution_t*) malloc(ga->popSize * sizeof(GaSolution_t));

    for (uint64_t i = 0; i < ga->popSize; i++) {
        ga->population[i] = GaSolution(ga->numberOfGenes, true);
    }

    ga->bestSolution = GaSolution(ga->numberOfGenes, true);
    ga->bestSolution.fitness = ga->fitnessEvaluationFn(ga->bestSolution.genes);
}

void RunGeneticAlgorithm(GeneticAlgorithm_t* const ga, const enum bool enablePrintOutput) {
	if (ga->block) {
		printf("\nGenetic algorithm not initialised properly. Will not run.");
		return;
	}
    
	CreatePopulation(ga);

	for (uint32_t i = 0; i < ga->numberOfGenerations; i++) {
		EvaluatePopFitness(ga);
		if (enablePrintOutput) printf("Generation %u completed. Best fitness = %f\n", i + 1, ga->bestSolution.fitness);
		TournamentSelection(ga);
	}
}

void EvaluatePopFitness(GeneticAlgorithm_t* const ga) {
    // this could be easily parrallised! TODO: multithread
    for (uint64_t i = 0; i < ga->popSize; i++) {
        ga->population[i].fitness = ga->fitnessEvaluationFn(
            ga->population[i].genes
        );

        if (ga->population[i].fitness > ga->bestSolution.fitness) {
            ga->bestSolution = CopyGaSolution(&(ga->population[i]));
        }
    }
}

void TournamentSelection(GeneticAlgorithm_t* const ga) {
    // this may be parallelised? TODO: multithread
    GaSolution_t* nextGen = (GaSolution_t*) malloc(ga->popSize * sizeof(GaSolution_t));    
    const GaSolution_t *mother, *father;

    for (uint32_t i = 0; i < ga->popSize; i++) {
        mother = GetBestContender(ga, NULL);
        father = GetBestContender(ga, mother);

        nextGen[i] = GenerateChild(ga, mother, father);
    }

    FreePopulation(ga);
    ga->population = nextGen;
    
}

const GaSolution_t* GetBestContender(const GeneticAlgorithm_t* const ga, const GaSolution_t* const dontSelect) {
    // Array of pointers to contenders
    const GaSolution_t** contenders = (const GaSolution_t**) malloc(ga->groupSize * sizeof(GaSolution_t*));
    const GaSolution_t* bestContender;
    enum bool isUniqueContender;

    // clear the list of contenders
    for (uint32_t i = 0; i < ga->groupSize; i++) {
        contenders[i] = NULL;
    }

    // get <groupSize> unique contenders
    for (uint32_t i = 0; i < ga->groupSize; i++) {
        isUniqueContender = false;
        GaSolution_t* randomSolution = NULL;
        while (isUniqueContender == false) {
            uint32_t randomSolutionIndex = GenerateRandomU32() % ga->popSize;
            randomSolution = &(ga->population[randomSolutionIndex]);

            // check that the contender isn't the same as in the "dontSelect" argument (i.e. already selected mother)
            isUniqueContender = !(randomSolution == dontSelect);

            // check that the contender isn't the same as any of the other contenders selected
            for (uint32_t j = 0; isUniqueContender == true && j < ga->groupSize; j++) {
                if (contenders[j] == randomSolution) isUniqueContender = false;
            }
        }
        if (randomSolution != NULL) contenders[i] = randomSolution;
    }

    // of the <groupSize> random contenders we just selected, find the best
    bestContender = contenders[0];
	for (uint32_t i = 0; i < ga->groupSize; i++) {
		if ((contenders[i]->fitness) > (bestContender->fitness)) {
			bestContender = contenders[i];
		}
	}

    free(contenders);
    return bestContender;
}

uint32_t GenerateRandomU32(void) {
    // existing randomisation functions weren't being random enough, not even <random>
	// so this function should be a little more random

	uint32_t u32RandomNumber = 0;

	for (uint64_t i = 0; i < 2; i++) {
		const uint16_t u16IntermediateRandom = rand() % 65536;
		const uint64_t u8Offset = i * 16;
		u32RandomNumber += (u16IntermediateRandom << u8Offset);
	}

	return u32RandomNumber;
}

byte GenerateRandomByte(void) {
	return rand() % 256;
}