#include "genalg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void FreeGeneticAlgorithm(GeneticAlgorithm_t* const ga) {
    free(ga->populationGenes);
    free(ga->populationFitness);
    free(ga->bestSolution.genes);
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
void GenerateChild(const GeneticAlgorithm_t* const ga, byte* child, const byte* const mother, const byte* const father) {
    uint32_t bitPosition;
    uint32_t crossMask = 0;
    uint32_t crossIndex = ga->numberOfGenes + 1;
    uint32_t mutateMask = 0;
    uint32_t mutateIndex = ga->numberOfGenes + 1;
    const byte bitsPerGene = 8;

    if (GenerateRandomByte() % 100 < ga->crossoverRate) {
        // crossover
        bitPosition = 1 + (GenerateRandomU32() % (bitsPerGene * ga->numberOfGenes - 1)); // N between 1 and bitsPerGene*NumberOfGenes - 1 (i.e. bit position between first and last bit of u32 array)
        crossIndex = bitPosition / bitsPerGene; // get the index of the byte array in which the bit must be
		crossMask = (1 << (bitPosition % bitsPerGene)) - 1; // get bitmask which will be all zeroes on 1 side and all 1s of the otherside, of the bit position
    }

    if (GenerateRandomByte() % 100 < ga->mutationRate) {
        //mutation
		bitPosition = 1 + (GenerateRandomU32() % (bitsPerGene * ga->numberOfGenes - 1)); // N between 1 and bitsPerGene*NumberOfGenes - 1 (i.e. bit position between first and last bit of u32 array)

		mutateIndex = bitPosition / bitsPerGene; // get the index of the byte array in which the bit must be
		mutateMask = (1 << (bitPosition % (bitsPerGene - 1))); // get bitmask which will be all zeroes on 1 side and all 1s of the otherside, of the bit position
    }

    for (uint32_t i = 0; i < ga->numberOfGenes; i++) {
		// crossover (aka where the mother's DNA meets the father's ... no, not like that)
		if (i < crossIndex)
			child[i] = mother[i];
		else if (i > crossIndex)
			child[i] = father[i];
		else
			child[i] = (mother[i] & crossMask) + (father[i] & ~crossMask);

		// mutation (aka flip a single bit)
		if (i == mutateIndex)
			child[i] = (child[i] & ~mutateMask) + (~(child[i]) & mutateMask);
	}
}

void CreatePopulation(GeneticAlgorithm_t* const ga) {
    const uint64_t popGenesMemSize = ga->popSize * ga->numberOfGenes * sizeof(byte);
    ga->populationGenes = (byte*) malloc(popGenesMemSize);
    ga->populationFitness = (double*) malloc(ga->popSize * sizeof(double));


    for (uint64_t i = 0; i < popGenesMemSize; i++) {
        ga->populationGenes[i] = GenerateRandomByte();
    }

    ga->bestSolution.genes = (byte*) malloc(ga->numberOfGenes * sizeof(byte));
    for (uint32_t i = 0; i < (ga->numberOfGenes * sizeof(byte)); i++) {
        ga->bestSolution.genes[i] = GenerateRandomByte();
    }

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
    for (uint32_t i = 0; i < ga->popSize; i++) {
        ga->populationFitness[i] = ga->fitnessEvaluationFn(
            &(ga->populationGenes[i * ga->numberOfGenes])
        );

        if (ga->populationFitness[i] > ga->bestSolution.fitness) {
            for (uint32_t j = 0; j < ga->numberOfGenes; j++) {
                ga->bestSolution.genes[j] = ga->populationGenes[i * ga->numberOfGenes + j];
                ga->bestSolution.fitness = ga->populationFitness[i];
            }
        }
    }
}

void TournamentSelection(GeneticAlgorithm_t* const ga) {
    // this may be parallelised? TODO: multithread
    byte* nextGenGenes = (byte*) malloc(ga->popSize * ga->numberOfGenes * sizeof(byte));    
    const byte *mother, *father;

    for (uint32_t i = 0; i < ga->popSize; i++) {
        mother = GetBestContenderGenes(ga, NULL);
        father = GetBestContenderGenes(ga, mother);

        GenerateChild(ga, &(nextGenGenes[i * ga->numberOfGenes]), mother, father);
    }

    free(ga->populationGenes);
    ga->populationGenes = nextGenGenes;
    
}

const byte* GetBestContenderGenes(const GeneticAlgorithm_t* const ga, const byte* const dontSelect) {
    uint32_t* contenders = (uint32_t*) malloc(ga->groupSize * sizeof(uint32_t));
    uint32_t bestContender;
    enum bool isUniqueContender;

    // get <groupSize> unique contenders
    for (uint32_t i = 0; i < ga->groupSize; i++) {
        isUniqueContender = false;
        uint32_t randomSolution;
        while (isUniqueContender == false) {
            randomSolution = GenerateRandomU32() % ga->popSize;
            // check that the contender isn't the same as in the "dontSelect" argument (i.e. already selected mother)
            const byte* const randomSolutionGeneAddr = &(ga->populationGenes[randomSolution * ga->numberOfGenes]);
            isUniqueContender = !(randomSolutionGeneAddr == dontSelect);

            // check that the contender isn't the same as any of the other contenders selected
            for (uint32_t j = 0; isUniqueContender == true && j < ga->groupSize; j++) {
                if (contenders[j] == randomSolution) isUniqueContender = false;
            }
        }
        contenders[i] = randomSolution;
    }

    // of the <groupSize> random contenders we just selected, find the best
    bestContender = contenders[0];
	for (uint32_t i = 0; i < ga->groupSize; i++) {
		if ((ga->populationFitness[contenders[i]]) > (ga->populationFitness[bestContender])) {
			bestContender = contenders[i];
		}
	}

    free(contenders);
    return &(ga->populationGenes[bestContender * ga->numberOfGenes]);
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