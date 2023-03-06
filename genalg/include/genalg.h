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
    double (*fitnessEvaluationFn)(const byte * const genes);
    uint32_t popSize;
    uint32_t mutationRate;
	uint32_t crossoverRate;
    uint16_t groupSize;
    uint16_t numberOfGenes;
    uint16_t numberOfGenerations;
    enum bool block;

    GaSolution_t* population;
    GaSolution_t bestSolution;
} GeneticAlgorithm_t;


void FreeGaSolution(GaSolution_t* const solution);
void MoveGaSolution(GaSolution_t* const from, GaSolution_t* const to);
GaSolution_t CopyGaSolution (const GaSolution_t* const from);
GENALG_IMEX GaSolution_t GaSolution(const uint16_t numberOfGenes, const enum bool isRandomGenes);

void FreePopulation(GeneticAlgorithm_t* const ga);
GENALG_IMEX void FreeGeneticAlgorithm(GeneticAlgorithm_t* const ga);
GENALG_IMEX void MoveGeneticAlgorithm(GeneticAlgorithm_t* const from, GeneticAlgorithm_t* const to);
GENALG_IMEX GeneticAlgorithm_t CopyGeneticAlgorithm(const GeneticAlgorithm_t* const from);
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

uint32_t GenerateRandomU32();
byte GenerateRandomByte();
GaSolution_t GenerateChild(const GeneticAlgorithm_t* const ga, const GaSolution_t* const mother, const GaSolution_t* const father);
void CreatePopulation(GeneticAlgorithm_t* const ga);

void EvaluatePopFitness(GeneticAlgorithm_t* const ga);
void TournamentSelection(GeneticAlgorithm_t* const ga);

const GaSolution_t* GetBestContender(const GeneticAlgorithm_t* const ga, const GaSolution_t* const dontSelect);
