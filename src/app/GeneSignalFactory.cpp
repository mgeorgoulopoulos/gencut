#include "GeneSignalFactory.h"

#include "GeneSignalList.h"
#include "GeneSignalMatrix.h"

#include "StringConstants.h"

#include <utils/Exception.h>

GeneSignalFactory::GeneSignalFactory() {}
GeneSignalFactory::~GeneSignalFactory() {}

std::vector<GeneSignalFactory::GeneSignalIdentifier>
GeneSignalFactory::availableSignalTypes() const {
	std::vector<GeneSignalIdentifier> result = {kmatrix, klist};
	return result;
}

GeneSignalPtr GeneSignalFactory::create(const GeneSignal::Options &options, AtomBox &geneRegistry) {
	GeneSignal *result = nullptr;

	if (options.signalType == kmatrix) {
		result = static_cast<GeneSignal *>(new GeneSignalMatrix(options, geneRegistry));
	}

	if (options.signalType == klist) {
		result = static_cast<GeneSignal *>(
			new GeneSignalList(options, geneRegistry));
	}

	if (result == nullptr) {
		throw(Exception({"Failed to create gene signal:", options.signalType}));
	}

	return std::shared_ptr<GeneSignal>(result);
}