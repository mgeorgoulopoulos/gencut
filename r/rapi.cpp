#ifndef SKIP_R_API

using namespace Rcpp;

// [[Rcpp::export]]
List gencutMatrix(NumericMatrix modelMatrix,
						   NumericMatrix signalMatrix,
						   std::string tailSelection) {

	CharacterVector clusters;
	NumericMatrix stats;

	try {
		AtomBox geneRegistry;

		GeneSignalPtr signal(
			new GeneSignalMatrix(signalMatrix, kaverage, geneRegistry));
		GenomeModel model(modelMatrix, geneRegistry);

		// Handshake model & signal by removing from the model nonexistent genes
		// in signal
		const GeneSet genesInSignal = signal->geneSet();
		model.removeGenesNotExistingIn(genesInSignal, geneRegistry);

		GenomeCutter::Options options;
		options.sphereRadius = 15.0;
		options.minimumGeneCount = 50;
		options.randomSphereCount = 10000;
		options.distributionExtractionIterations = 10000;
		options.pAdjThreshold = 0.01;
		options.overlapThreshold = 0.05;
		options.tailSelection =
			GenomeCutter::Options::tailSelectionFromString(tailSelection);

		GenomeCutter cutter(model, signal, options);
		GenomeCutter::Result result = cutter.cut();

		// Gen uncustered genes
		GeneSet unclusteredGenes = model.geneSet();
		for (const GeneSet &cluster : result.clusters) {
			for (GeneId gene : cluster) {
				unclusteredGenes.erase(gene);
			}
		}

		// Write output
		for (int i = 0; i < (int)result.clusters.size(); i++) {
			const std::string groupLetter = std::string(1, 'A' + i);
			for (const GeneId gene : result.clusters[i]) {
				const std::string geneName = geneRegistry.name(gene);
				clusters.push_back(groupLetter, geneName);
			}
		}
		// And the unclustered ones, with <null> cluster name
		for (const GeneId gene : unclusteredGenes) {
			const std::string geneName = geneRegistry.name(gene);
			clusters.push_back("none", geneName);
		}

		// Output statistics
		stats = NumericMatrix(result.sampleStatistics.size(), 4);
		colnames(stats) = CharacterVector::create("Metric", "Shuffled", "pValue", "pAdj");
		for (int i = 0; i < (int)result.sampleStatistics.size(); i++) {
			const GenomeCutter::SampleStats &s = result.sampleStatistics[i];
			stats(i, 0) = s.metric;
			stats(i, 1) = s.shuffled;
			stats(i, 2) = s.pValue;
			stats(i, 3) = s.adjustedPValue;
		}

	} catch (Exception e) {
		e.print();
		throw "gencut error";
	}

	return List::create(Named("clusters") = clusters, Named("stats") = stats);
}

/*** R

*/

#endif // ndef SKIP_R_API
