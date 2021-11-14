#ifndef __GENOME_CUTTER_H__
#define __GENOME_CUTTER_H__

#include "GeneCollections.h"
#include "GeneSignal.h"

#include <vector>


class GenomeModel;

class GenomeCutter {
  public:

	struct Options {
		// Radius of the sampling sphere, in model space.
		double sphereRadius = 15.0;

		// Threshold to disregard mostly empty spheres
		int minimumGeneCount = 50;

		// How many successful samples to consider
		int randomSphereCount = 10000;

		// How many samples to take in order to assess the distribution of the
		// metric. The more, the greater resolution we get in p-values.
		int distributionExtractionIterations = 10000;

		// Adjusted p-value threshold for considering a sphere-sample
		// "interesting".
		double pAdjThreshold = 0.01;

		// Used during hierarchical clustering of the "interesting"
		// sphere-samples. We consider two clusters as distinct when they
		// overlap less than this.
		double overlapThreshold = 0.05;

		// Which tail (or both) of the distribution do we consider "interesting"?
		enum class TailSelection {Low, High, Both};
		TailSelection tailSelection = TailSelection::Low;

		void print() const;
		void printHelp() const;
	};

	GenomeCutter(const GenomeModel &model, const GeneSignalPtr signal,
				 const Options &options);
	virtual ~GenomeCutter();

	// Statistics of a sphere-sample.
	struct SampleStats {
		double metric = 0.0;
		double pValue = 0.0;
		double adjustedPValue = 0.0;
	};

	struct Result {
		std::vector<GeneSet> clusters;
		std::vector<SampleStats> sampleStatistics;
	};

	// Produces segments of the gene model using signal
	Result cut() const;

  private:
	struct Private;
	Private *d = nullptr;
};

#endif // __GENOME_CUTTER_H__