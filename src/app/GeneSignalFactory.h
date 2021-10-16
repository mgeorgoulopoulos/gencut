#ifndef __GENE_SIGNAL_FACTORY_H__
#define __GENE_SIGNAL_FACTORY_H__

#include "GeneSignal.h"

class GeneSignalFactory {
  public:
	GeneSignalFactory();
	virtual ~GeneSignalFactory();

	using GeneSignalIdentifier = std::string;

	std::vector<GeneSignalIdentifier> availableSignalTypes() const;

	GeneSignalPtr create(const GeneSignal::Options &options,
						 AtomBox &geneRegistry);
};

#endif // __GENE_SIGNAL_FACTORY_H__