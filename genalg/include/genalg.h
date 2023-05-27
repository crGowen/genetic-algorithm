#pragma once

#ifdef _WIN32
	#ifdef GENALG_EXPORT
	#define GENALG_IMEX __declspec(dllexport)
	#else
	#define GENALG_IMEX __declspec(dllimport)
	#endif
#else
	#define GENALG_IMEX
#endif

#include <stdint.h>

enum bool {false, true};
typedef unsigned char byte;

GENALG_IMEX typedef struct {
    double fitness;
    byte* genes;
} BestSolution_t;

GENALG_IMEX typedef struct {
    BestSolution_t bestSolution;
    double (*fitnessEvaluationFn)(const byte * const genes);
    double* populationFitness;
    byte* populationGenes;
    uint32_t popSize;
    enum bool block;
    uint16_t groupSize;
    uint16_t numberOfGenes;
    uint16_t numberOfGenerations;
    uint8_t mutationRate;
	uint8_t crossoverRate;
} GeneticAlgorithm_t;

GENALG_IMEX void FreeGeneticAlgorithm(GeneticAlgorithm_t* const ga);
GENALG_IMEX GeneticAlgorithm_t GeneticAlgorithm(
    const uint32_t popSize,
    const uint8_t mutationRate,
    const uint8_t crossoverRate,
    const uint16_t groupSize,
    const uint16_t numberOfGenes,
    const uint16_t numberOfGenerations,
    double (* const fitnessEvaluationFn)(const byte* const genes)
);

GENALG_IMEX void RunGeneticAlgorithm(GeneticAlgorithm_t* const ga, const uint8_t numThreads, const enum bool enablePrintOutput);