#ifndef __GENE_SIGNAL_MATRIX_H__
#define __GENE_SIGNAL_MATRIX_H__

#include "GeneSignal.h"

// GeneSignal implementation based on a matrix of gene->gene scores.
class GeneSignalMatrix : public GeneSignal {
  public:
	GeneSignalMatrix(const Options &options, AtomBox &geneRegistry);
	virtual ~GeneSignalMatrix();

	virtual GeneSet geneSet() const override;

	double metric(const GeneList &genes) const override;

	double cell(GeneId a, GeneId b) const;

private:
	struct Private;
  Private *d = nullptr;

};

#endif // __GENE_SIGNAL_MATRIX_H__