#include "GeneSignal.h"

void GeneSignal::Options::print() const {
	printf("\tSignal type: %s\n", signalType.c_str());
	printf("\tMetric: %s\n", metric.c_str());
}

void GeneSignal::Options::printHelp() const {
	printf("\t--signal filename\t\tLoad gene signal from CSV file 'filename'.\n");
	printf("\t--signal-type type\t\tSpecify type of signal. Currently only 'matrix' is valid.\n");
	printf("\t--metric function\t\tSpecify function to be applied to sets of genes, ex: 'average'.\n");
}

GeneSignal::GeneSignal(const Options &options, AtomBox &geneRegistry)
	: options(options), geneRegistry(geneRegistry) {}
GeneSignal::~GeneSignal() {}