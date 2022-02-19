#ifndef __GENE_SIGNAL_LIST_H__
#define __GENE_SIGNAL_LIST_H__

#include "GeneSignal.h"

// GeneSignal implementation based on a list of genes and properties. Properties may be numeric or factors.
class GeneSignalList : public GeneSignal {
  public:
	GeneSignalList(const Options &options, AtomBox &geneRegistry);
	virtual ~GeneSignalList();

	virtual GeneSet geneSet() const override;

	double metric(const GeneList &genes) const override;

private:
	struct Private;
  Private *d = nullptr;

};

#endif // __GENE_SIGNAL_LIST_H__