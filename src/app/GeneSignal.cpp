#include "GeneSignal.h"

void GeneSignal::Options::print() const {
	printf("\tSignal type: %s\n", signalType.c_str());
	printf("\tMetric: %s\n", metric.c_str());
}

GeneSignal::GeneSignal(const Options &options, AtomBox &geneRegistry)
	: options(options), geneRegistry(geneRegistry) {}
GeneSignal::~GeneSignal() {}