#ifndef __GENE_SIGNAL_H__
#define __GENE_SIGNAL_H__

#include "GeneCollections.h"

#include <memory>

// Abstract class. Calculates metric for sets of genes.
class GeneSignal {
  public:
	struct Options {
		std::string dataFilename;
		std::string signalType;
		std::string metric;

		void print() const;
		void printHelp() const;
	};
	GeneSignal(const Options &options, AtomBox &geneRegistry);
	virtual ~GeneSignal();

	virtual GeneSet geneSet() const = 0;

	virtual double metric(const GeneList &genes) const = 0;

  private:
	Options options;
	AtomBox &geneRegistry;
};

using GeneSignalPtr = std::shared_ptr<GeneSignal>;

#endif // __GENE_SIGNAL_H__