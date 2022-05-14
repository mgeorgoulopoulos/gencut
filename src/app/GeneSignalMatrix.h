#ifndef __GENE_SIGNAL_MATRIX_H__
#define __GENE_SIGNAL_MATRIX_H__

#include "GeneSignal.h"

// GeneSignal implementation based on a matrix of gene->gene scores.
class GeneSignalMatrix : public GeneSignal {
  public:
#ifndef GENCUT_R_PACKAGE
	GeneSignalMatrix(const Options &options, AtomBox &geneRegistry);
#endif // ndef GENCUT_R_PACKAGE

#ifndef SKIP_R_API
	// Allow to create a matrix from Rcpp::NumericMatrix
	GeneSignalMatrix(Rcpp::NumericMatrix m, const std::string &metric, AtomBox &geneRegistry);
#endif // ndef SKIP_R_API

	virtual ~GeneSignalMatrix();

	virtual GeneSet geneSet() const override;

	double metric(const GeneList &genes) const override;

	double cell(GeneId a, GeneId b) const;

private:
	struct Private;
  Private *d = nullptr;

};

#endif // __GENE_SIGNAL_MATRIX_H__