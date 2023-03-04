#include "genalg.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// GaSolution constructors/destructors
void FreeGaSolution(GaSolution_t* const solution) {
    if (solution->genes != NULL) {
        free(solution->genes);
    }
}

void MoveGaSolution(GaSolution_t* const from, GaSolution_t* const to) {
    FreeGaSolution(to);

    to->fitness = from->fitness;
    to->numberOfGenes = from->numberOfGenes;
    to->genes = from->genes;
    from->genes = NULL;
}

GaSolution_t CopyGaSolution (const GaSolution_t* const from) {
    GaSolution_t to;

    to.fitness = from->fitness;
    to.numberOfGenes = from->numberOfGenes;

    to.genes = (uint32_t*) malloc(to.numberOfGenes * sizeof(uint32_t));
    for (uint64_t i = 0; i < to.numberOfGenes; i++) {
            to.genes[i] = from->genes[i];
    }
    
    return to;
}

GaSolution_t GaSolution(const uint16_t numberOfGenes, const enum bool isRandomGenes) {
    GaSolution_t solution;
    solution.numberOfGenes = numberOfGenes;
    solution.genes = (uint32_t*) malloc(numberOfGenes * sizeof(uint32_t));

    if (isRandomGenes == true) {
        for (uint64_t i = 0; i < numberOfGenes; i++) {
            solution.genes[i] = GenerateRandomU32();
        }
    }

    solution.fitness = 0.0;

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

void MoveGeneticAlgorithm(GeneticAlgorithm_t* const from, GeneticAlgorithm_t* const to) {
    FreeGeneticAlgorithm(to);

    to->popSize                 = from->popSize;
    to->mutationRate            = from->mutationRate;
    to->crossoverRate           = from->crossoverRate;
    to->groupSize               = from->groupSize;
    to->numberOfGenes           = from->numberOfGenes;
    to->numberOfGenerations     = from->numberOfGenerations;
    to->fitnessEvaluationFn     = from->fitnessEvaluationFn;

    to->population = from->population;
    from->population = NULL;

    to->block = from->block;

    MoveGaSolution(&(from->bestSolution), &(to->bestSolution));
}

GeneticAlgorithm_t CopyGeneticAlgorithm(const GeneticAlgorithm_t* const from) {
    GeneticAlgorithm_t to;

    to.popSize                 = from->popSize;
    to.mutationRate            = from->mutationRate;
    to.crossoverRate           = from->crossoverRate;
    to.groupSize               = from->groupSize;
    to.numberOfGenes           = from->numberOfGenes;
    to.numberOfGenerations     = from->numberOfGenerations;
    to.fitnessEvaluationFn     = from->fitnessEvaluationFn;

    to.population = (GaSolution_t*) malloc(to.popSize * sizeof(GaSolution_t));
    for (uint64_t i = 0; i < to.popSize; i++) {
            to.population[i] = CopyGaSolution(&(from->population[i]));
    }

    to.block = from->block;

    to.bestSolution = CopyGaSolution(&(from->bestSolution));

    return to;
}

GeneticAlgorithm_t GeneticAlgorithm(
    const uint32_t popSize,
    const uint32_t mutationRate,
    const uint32_t crossoverRate,
    const uint16_t groupSize,
    const uint16_t numberOfGenes,
    const uint16_t numberOfGenerations,
    double (* const fitnessEvaluationFn)(const uint32_t* const genes)
) {
    GeneticAlgorithm_t ga;
    ga.popSize              = popSize;
    ga.mutationRate         = mutationRate;
    ga.crossoverRate        = crossoverRate;
    ga.groupSize            = groupSize;
    ga.numberOfGenes        = numberOfGenes;
    ga.numberOfGenerations  = numberOfGenerations;
    ga.fitnessEvaluationFn  = fitnessEvaluationFn;

    ga.population = NULL;

    ga.block = false;

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

    if (GenerateRandomU32() % 100000 < ga->crossoverRate) {
        // crossover
        bitPosition = 1 + (GenerateRandomU32() % (32 * ga->numberOfGenes - 1)); // N between 1 and 32*NumberOfGenes - 1 (i.e. bit position between first and last bit of u32 array)
        crossIndex = bitPosition / 32; // get the index of the u32 array in which the bit must be
		crossMask = (1 << (bitPosition % 32)) - 1; // get bitmask which will be all zeroes on 1 side and all 1s of the otherside, of the bit position
    }

    if (GenerateRandomU32() % 100000 < ga->mutationRate) {
        //mutation
		bitPosition = 1 + (GenerateRandomU32() % (32 * ga->numberOfGenes - 1)); // N between 1 and 32*NumberOfGenes - 1 (i.e. bit position between first and last bit of u32 array)

		mutateIndex = bitPosition / 32; // get the index of the u32 array in which the bit must be
		mutateMask = (1 << (bitPosition % (32 - 1))); // get bitmask which will be all zeroes on 1 side and all 1s of the otherside, of the bit position
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

    GaSolution_t t = GenerateChild(ga, &(ga->population[1]), &(ga->population[2]));

    MoveGaSolution(
        &t,
        &(ga->population[0])
        );

	for (uint32_t i = 0; i < ga->numberOfGenerations; i++) {
		EvaluatePopFitness(ga);
		if (enablePrintOutput) printf("Generation %u completed. Best fitness = %f\n", i + 1, ga->bestSolution.fitness);
		TournamentSelection(ga);
	}
}

void EvaluatePopFitness(GeneticAlgorithm_t* const ga) {
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

uint32_t GenerateRandomU32() {
    // existing randomisation functions weren't being random enough, not even <random>
	// so this function should be a little more random

	uint32_t u32RandomNumber = 0;

	for (uint64_t i = 0; i < 2; i++) {
		const uint16_t u16IntermediateRandom = rand() % 65536;
		const uint8_t u8Offset = i * 16;
		u32RandomNumber += (u16IntermediateRandom << u8Offset);
	}

	return u32RandomNumber;
}