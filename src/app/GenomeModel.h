#ifndef __GENOME_MODEL_H__
#define __GENOME_MODEL_H__

#include "GeneCollections.h"

#include <utils/AtomBox.h>
#include <utils/Vec3D.h>

#include <vector>

using ChromosomeId = Atom;

// A continuous segment of a chromosome
struct Locus {
	ChromosomeId chromosomeId = 0;

	// Start and end nucleotide positions (inclusive).
	int startPosition = 0;
	int endPosition = 0;
};

// Represents a 3D model of a genome.
class GenomeModel {
  public:
	struct Gene {
		GeneId geneId = 0;
		Locus locus;
		Vec3D position;
	};

	struct Options {
		std::string filename;

		void printHelp() const;
	};

	GenomeModel();
	
#ifndef GENCUT_R_PACKAGE
	GenomeModel(const Options &options, AtomBox &geneRegistry);
#endif // ndef GENCUT_R_PACKAGE

#ifndef SKIP_R_API
	GenomeModel(Rcpp::NumericMatrix m, AtomBox &geneRegistry);
#endif // ndef SKIP_R_API


	GeneSet geneSet() const;

	GeneList geneList() const;

	void removeGenes(const GeneSet &genesToRemove);

	GeneList sphereSample(const Vec3D &center, double radius) const;

	Vec3D minBound() const;
	Vec3D maxBound() const;

	void removeGenesNotExistingIn(const GeneSet &genes, AtomBox &geneRegistry);

  private:
	std::vector<Gene> genes;
};

#endif // __GENOME_MODEL_H__