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
    uint16_t numberOfGenes;
} GaSolution_t;

GENALG_IMEX typedef struct {
    GaSolution_t bestSolution;
    double (*fitnessEvaluationFn)(const byte * const genes);
    GaSolution_t* population;
    uint32_t popSize;
    uint32_t mutationRate;
	uint32_t crossoverRate;
    enum bool block;
    uint16_t groupSize;
    uint16_t numberOfGenes;
    uint16_t numberOfGenerations;
} GeneticAlgorithm_t;


void FreeGaSolution(GaSolution_t* const solution);
GaSolution_t CopyGaSolution (const GaSolution_t* const from);
GENALG_IMEX GaSolution_t GaSolution(const uint16_t numberOfGenes, const enum bool isRandomGenes);

void FreePopulation(GeneticAlgorithm_t* const ga);
GENALG_IMEX void FreeGeneticAlgorithm(GeneticAlgorithm_t* const ga);
GENALG_IMEX GeneticAlgorithm_t GeneticAlgorithm(
    const uint32_t popSize,
    const uint32_t mutationRate,
    const uint32_t crossoverRate,
    const uint16_t groupSize,
    const uint16_t numberOfGenes,
    const uint16_t numberOfGenerations,
    double (* const fitnessEvaluationFn)(const byte* const genes)
);

GENALG_IMEX void RunGeneticAlgorithm(GeneticAlgorithm_t* const ga, const enum bool enablePrintOutput);

uint32_t GenerateRandomU32(void);
byte GenerateRandomByte(void);
GaSolution_t GenerateChild(const GeneticAlgorithm_t* const ga, const GaSolution_t* const mother, const GaSolution_t* const father);
void CreatePopulation(GeneticAlgorithm_t* const ga);

void EvaluatePopFitness(GeneticAlgorithm_t* const ga);
void TournamentSelection(GeneticAlgorithm_t* const ga);

const GaSolution_t* GetBestContender(const GeneticAlgorithm_t* const ga, const GaSolution_t* const dontSelect);
