#include "GenomeCutter.h"

#include "GeneSignal.h"
#include "GenomeModel.h"

#include <utils/Exception.h>

#include <omp.h>

#include <iterator>
#include <random>

void GenomeCutter::Options::print() const {
	printf("\tSphere radius: %f\n", sphereRadius);
	printf("\tMinimum gene count: %d\n", minimumGeneCount);
	printf("\tRandom sample count: %d\n", randomSphereCount);
	printf("\tDistribution extraction iterations: %d\n",
		   distributionExtractionIterations);
	printf("\tAdjusted p-value threshold: %f\n", pAdjThreshold);
	printf("\tOverlap threshold: %f\n", overlapThreshold);

	printf("\tTail selection: ");
	switch (tailSelection) {
	case TailSelection::Low:
		printf("low");
		break;
	case TailSelection::High:
		printf("high");
		break;
	case TailSelection::Both:
		printf("both");
		break;
	}
	printf("\n");
}

void GenomeCutter::Options::printHelp() const {
	printf("\t--radius r\t\tSpecify radius of the sampling sphere, in model "
		   "space. Default is r=15.0\n");
	printf("\t--min-genes N\t\tAn acceptable sphere-sample must have at least "
		   "N genes. Default is 50.\n");
	printf("\t--samples N\t\tGenerate N random spheres. Default is 10^5.\n");
	printf(
		"\t--p-samples N\t\tCalculate the metric over N random sets of genes "
		"(not spatially associated) to assess signal distribution. Default is "
		"10^5 samples which yields a p-value resolution of 10^-5.\n");
	printf("\t--padj-threshold p\t\tThreshold of adjusted p-value to be "
		   "considered as significant. Default is p=0.01.\n");
	printf("\t--overlap-threshold t\t\tOverlap raio considered as distinct. "
		   "Default is 0.05, which means two gene sets A and B are considered "
		   "distinct when the share less than 5%% of their genes.\n");
	printf("\t--tail low|high|both\t\tWhich end of the distribution will be "
		   "used to calculate a p-value? 'both' specifies a double-tailed "
		   "test. 'low' and 'high' specify single-tailed.\n");
}

namespace {
struct WorkUnit {
	GeneList genes;
	double metric = 0.0;
	double pValue = 1.0;

	// P-value rank and Benjamini-Hochberg-adjusted version of it.
	int rank = 0;
	double adjustedPValue = 1.0;
};

// Adjusts p-values and reorders all work units.
void benjamini(std::vector<WorkUnit> &workUnits) {
	// Sort p-values from smaller to larger
	auto pValueCompare = [](const WorkUnit &a, const WorkUnit &b) {
		return a.pValue < b.pValue;
	};
	std::sort(workUnits.begin(), workUnits.end(), pValueCompare);

	// Apply Benjamini-Hochberg correction
	double previousPValue = workUnits.back().pValue;
	for (int i = (int)workUnits.size() - 1; i >= 0; i--) {
		WorkUnit &workUnit = workUnits[i];
		workUnit.rank = i + 1;
		workUnit.adjustedPValue =
			std::min(previousPValue, workUnit.pValue * workUnits.size() /
										 (double)workUnit.rank);
		previousPValue = workUnit.adjustedPValue;
	}
}

// Given a list of work units, it combines the overlapping spheres into clusters
std::vector<GeneSet>
clusterByGeneOverlap(const std::vector<WorkUnit> &workUnits,
					 double overlapThreshold, double *maximumOverlapRatio) {
	std::vector<GeneSet> result;
	for (const WorkUnit &workUnit : workUnits) {
		GeneSet geneSet;
		for (const GeneId gene : workUnit.genes) {
			geneSet.insert(gene);
		}
		result.push_back(geneSet);
	}

	// Now do hierarchical clustering.
	while (true) {
		// Find max-overlap pair of clusters
		*maximumOverlapRatio = 0.0;
		using ClusterPair = std::pair<int, int>;
		ClusterPair bestClusterPair(0, 1);
		GeneSet mergedCluster;
		for (int i = 0; i < (int)result.size(); i++) {
			const GeneSet &clusterA = result[i];
			for (int j = i + 1; j < (int)result.size(); j++) {
				const GeneSet &clusterB = result[j];
				const int minSize =
					std::min((int)clusterA.size(), (int)clusterB.size());
				GeneSet intersection;
				std::set_intersection(
					clusterA.begin(), clusterA.end(), clusterB.begin(),
					clusterB.end(),
					std::inserter(intersection, intersection.begin()));

				const int overlapCount = (int)intersection.size();
				const double overlapRatio =
					(double)overlapCount / (double)minSize;
				if (overlapRatio > *maximumOverlapRatio) {
					*maximumOverlapRatio = overlapRatio;

					mergedCluster.clear();
					std::set_union(
						clusterA.begin(), clusterA.end(), clusterB.begin(),
						clusterB.end(),
						std::inserter(mergedCluster, mergedCluster.begin()));
					bestClusterPair = ClusterPair(i, j);
				}
			}
		} // end for (all pairs)

		// When we reach the point where the best-overlapping clusters may be
		// considered distinct, we are done.
		if (*maximumOverlapRatio < overlapThreshold) {
			break;
		}

		// Not done yet: apply the merge and continue.
		result[bestClusterPair.second] = result.back();
		result.pop_back();
		result[bestClusterPair.first] =
			result.back(); // we can do this always: first is less than second
						   // by design.
		result.pop_back();
		result.push_back(mergedCluster);
	} // end hierarchical clustering

	return result;
}

} // end anonymous namespace

struct GenomeCutter::Private {
	Private(const GenomeModel &model, const GeneSignalPtr signal,
			const Options &options)
		: model(model), signal(signal), options(options),
		  randomGenerator(std::random_device()()) {

		genes = model.geneList();
		geneIndexDistribution =
			std::uniform_int_distribution<int>(0, (int)genes.size() - 1);

		// Get bounds of the model
		modelMinBound = model.minBound();
		modelMaxBound = model.maxBound();

		// Initialize distributions based on bounds
		xDistribution = std::uniform_real_distribution<double>(modelMinBound.x,
															   modelMaxBound.x);
		yDistribution = std::uniform_real_distribution<double>(modelMinBound.y,
															   modelMaxBound.y);
		zDistribution = std::uniform_real_distribution<double>(modelMinBound.z,
															   modelMaxBound.z);
	}

	double
	calculatePValue(double metric,
					const std::vector<double> &metricInRandomSamples) const;

	// Create a work unit with a set of randomly picked genes
	WorkUnit createRandomWorkUnit(int geneCount) const;

	// Create a work unit with genes sampled from a random sphere
	WorkUnit createSphereWorkUnit() const;

	const GenomeModel &model;
	const GeneSignalPtr signal;
	const Options options;

	GeneList genes;

	// For sampling random genes
	mutable std::default_random_engine randomGenerator;
	mutable std::uniform_int_distribution<int> geneIndexDistribution;

	// For sampling random-sphere genes
	Vec3D modelMinBound;
	Vec3D modelMaxBound;
	mutable std::uniform_real_distribution<double> xDistribution;
	mutable std::uniform_real_distribution<double> yDistribution;
	mutable std::uniform_real_distribution<double> zDistribution;
};

GenomeCutter::GenomeCutter(const GenomeModel &model, const GeneSignalPtr signal,
						   const Options &options)
	: d(new Private(model, signal, options)) {}

GenomeCutter::~GenomeCutter() { delete d; }

std::vector<GeneSet> GenomeCutter::cut() const {
	// Create n sphere samples
	printf("Sampling %d spheres... ", d->options.randomSphereCount);
	std::vector<WorkUnit> workUnits;
	std::set<int> setOfSphereGeneCounts;
	workUnits.reserve(d->options.randomSphereCount);
	for (int i = 0; i < d->options.randomSphereCount; i++) {
		WorkUnit workUnit = d->createSphereWorkUnit();
		workUnits.push_back(workUnit);

		// Keep track of this. It will come in handy later
		setOfSphereGeneCounts.insert((int)workUnit.genes.size());
	}
	printf("Done\n");

	// Parallelize the following loop.
	omp_lock_t writeLock;
	omp_init_lock(&writeLock);

	// Convert to list for OpenMP to consume
	std::vector<int> listOfGeneCounts;
	for (const int geneCount : setOfSphereGeneCounts)
		listOfGeneCounts.push_back(geneCount);

	// For each possible group size (ex: 50 genes, 118 genes, etc...), keep a
	// list of metric values coming from this many random gene collections of
	// that size.
	printf(
		"Assessing distribution of metric for %d possible gene set sizes ... ",
		setOfSphereGeneCounts.size());
	std::map<int, std::vector<double>> geneCountToDistribution;
#pragma omp parallel for
	for (int i = 0; i < (int)listOfGeneCounts.size(); i++) {
		const int geneCount = listOfGeneCounts[i];
		// Assess distribution of sets of 'geneCount' genes.
		std::vector<double> distribution;
		distribution.reserve(d->options.distributionExtractionIterations);
		for (int j = 0; j < d->options.distributionExtractionIterations; j++) {
			WorkUnit randomUnit = d->createRandomWorkUnit(geneCount);
			const double metric = d->signal->metric(randomUnit.genes);
			distribution.push_back(metric);
		}

		// Write shared state one thread at a time
		omp_set_lock(&writeLock);
		geneCountToDistribution[geneCount] = distribution;
		omp_unset_lock(&writeLock);
	}
	printf("Done\n");

	omp_destroy_lock(&writeLock);

	// Calculate metric and p-value for each sphere
	printf("Calculating metric and p-value for %d spheres ... ",
		   workUnits.size());
	for (WorkUnit &workUnit : workUnits) {
		workUnit.metric = d->signal->metric(workUnit.genes);
		workUnit.pValue = d->calculatePValue(
			workUnit.metric,
			geneCountToDistribution[(int)workUnit.genes.size()]);
	}
	printf("Done\n");

	// Adjust p-values
	printf("Adjusting p-values using Benjamini-Hochberg method... ");
	benjamini(workUnits);
	printf("Done.\n");

	// Filter by adjusted p-value
	std::vector<WorkUnit> tmp;
	std::copy_if(workUnits.begin(), workUnits.end(), std::back_inserter(tmp),
				 [&](const WorkUnit &x) {
					 return x.adjustedPValue <= d->options.pAdjThreshold;
				 });
	workUnits = tmp;
	printf("%d significant p-values (below %.05f)\n", workUnits.size(),
		   d->options.pAdjThreshold);

	// Get the significant genes
	GeneSet significantGenes;
	for (const WorkUnit &workUnit : workUnits) {
		for (const GeneId geneId : workUnit.genes) {
			significantGenes.insert(geneId);
		}
	}
	printf("%d significant genes.\n", significantGenes.size());
	;
	if (workUnits.empty()) {
		printf("No significant samples found");
	}

	// Convert work units to clusters of genes
	printf("\n");
	printf("Hierarchical clustering. Using threshold of %.02f%% overlap ratio "
		   "to consider clusters as distinct ... ",
		   d->options.overlapThreshold * 100.0);
	double maximumOverlapRatio = 0.0;
	std::vector<GeneSet> result = clusterByGeneOverlap(
		workUnits, d->options.overlapThreshold, &maximumOverlapRatio);
	printf("Done.\n");
	printf("Stopping clustering with %d clusters, %.02f%% maximum gene "
		   "overlap.\n",
		   result.size(), maximumOverlapRatio * 100.0);

	// The returned clusters are fairly defined. Sometimes we get a different
	// number of clusters because this is a random sampling, but mostly we get
	// the same. Ordering the clusters by size.
	std::sort(
		result.begin(), result.end(),
		[](const GeneSet &a, const GeneSet &b) { return a.size() < b.size(); });

	// Cleanup duplicate genes - let the smaller cluster retain all its genes
	GeneSet genesUsed;
	for (GeneSet &cluster : result) {
		GeneSet tmp;
		std::set_difference(cluster.begin(), cluster.end(), genesUsed.begin(),
							genesUsed.end(), std::inserter(tmp, tmp.begin()));
		cluster = tmp;

		tmp.clear();
		std::set_union(genesUsed.begin(), genesUsed.end(), cluster.begin(),
					   cluster.end(), std::inserter(tmp, tmp.begin()));
		genesUsed = tmp;
	}

	// Report clusters
	for (int i = 0; i < result.size(); i++) {
		printf("\tCluster %d: %d genes\n", i + 1, result[i].size());
	}

	// Done
	return result;
}

double GenomeCutter::Private::calculatePValue(
	double metric, const std::vector<double> &metricInRandomSamples) const {

	if (metricInRandomSamples.empty()) {
		throw(Exception("Empty random sample set."));
	}

	// See how many times random is more extreme than metric.
	int randomLowerCount = 0;
	int randomHigherCount = 0;
	for (const double metricInRandom : metricInRandomSamples) {
		if (metricInRandom <= metric)
			randomLowerCount++;
		if (metricInRandom >= metric)
			randomHigherCount++;
	}

	const double pValueLow =
		(double)randomLowerCount / (double)metricInRandomSamples.size();
	const double pValueHigh =
		(double)randomHigherCount / (double)metricInRandomSamples.size();

	double pValue = 1.0;
	switch (options.tailSelection) {
	case Options::TailSelection::Low:
		pValue = pValueLow;
		break;
	case Options::TailSelection::High:
		pValue = pValueHigh;
		break;
	case Options::TailSelection::Both:
		pValue = 2.0 * std::min(pValueHigh, pValueLow);
		pValue = std::max(1.0, pValue);
		break;
	}

	// Give benefit of the doubt to chance: replace zero pValues with the
	// smallest we can safely say
	pValue = std::max(pValue, 1.0 / (double)metricInRandomSamples.size());

	return pValue;
}

WorkUnit GenomeCutter::Private::createRandomWorkUnit(int geneCount) const {
	WorkUnit result;

	if (geneCount <= 0)
		return result;

	if (genes.empty()) {
		throw(Exception("Cannot sample empty gene set"));
	}

	for (int i = 0; i < geneCount; i++) {
		result.genes.push_back(genes[geneIndexDistribution(randomGenerator)]);
	}

	return result;
}

WorkUnit GenomeCutter::Private::createSphereWorkUnit() const {
	WorkUnit result;

	while ((int)result.genes.size() < options.minimumGeneCount) {
		// Pick a random center for the sphere
		Vec3D center;
		center.x = xDistribution(randomGenerator);
		center.y = yDistribution(randomGenerator);
		center.z = zDistribution(randomGenerator);

		result.genes = model.sphereSample(center, options.sphereRadius);
	}

	return result;
}