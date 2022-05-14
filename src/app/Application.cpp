#include "Application.h"

#include "GeneSignalFactory.h"
#include "GenomeCutter.h"
#include "GenomeModel.h"

#include <algorithm>
#include <chrono>
#include <iterator>

void Application::Arguments::printHelp() const {
	printf("gencut [--help] [--binaryArguments values]\n");
	printf("\t--help\t\tDisplay this help message and exist\n");
	printf("\n");
	printf("General options\n");
	printf("\t--settings filename\t\tLoad arguments from 'filename'. This is a "
		   "CSV file containing two columns: 'key' and 'value'. 'key' column "
		   "expects the same named arguments used in command line. For "
		   "example, using --radius 123.0 in commad line is equivalent to "
		   "loading a settings file with two rows: header will contain "
		   "'key,value' and the row below will contain 'radius,123.0'. Please "
		   "omit the '--' prefix in the CSV.\n");
	printf("\n");
	printf("Genome 3D model options:\n");
	modelOptions.printHelp();
	printf("\n");
	printf("Gene signal options:\n");
	signalOptions.printHelp();
	printf("\n");
	printf("Processing options:\n");
	cutterOptions.printHelp();
	printf("\n");
}

void Application::Arguments::print() const {
	printf("\tGenome 3D model: %s\n", modelOptions.filename.c_str());
	printf("\tGene signal: %s\n", signalOptions.dataFilename.c_str());
	printf("\tOutput file: %s\n", outputFilename.c_str());
	printf("\n");
	printf("Signal options:\n");
	signalOptions.print();
	printf("Segmentation options:\n");
	cutterOptions.print();
}

struct Application::Private {
  	Private(const Application::Arguments &arguments)
		: arguments(arguments) {}

	void initializeObjects();

	Application::Arguments arguments;

	AtomBox geneRegistry;
	GenomeModel model;
	GeneSignalPtr signal;
};

Application::Application(const Arguments &arguments)
	: d(new Private(arguments)) {}

Application::~Application() { delete d; }

void Application::execute() {
	if (d->arguments.displayHelp) {
		d->arguments.printHelp();
		return;
	}

	// Start clock
	std::chrono::steady_clock::time_point startTime =
		std::chrono::steady_clock::now();

	// Print options
	printf("Configuration:\n");
	d->arguments.print();
	printf("---------------------------\n");

	d->initializeObjects();

	// Finally ... cut!
	GenomeCutter cutter(d->model, d->signal, d->arguments.cutterOptions);
	GenomeCutter::Result result = cutter.cut();

	// Get set of genes not belonging to any cluster, to accompany the rest
	GeneSet unclusteredGenes = d->model.geneSet();
	for (const GeneSet &cluster : result.clusters) {
		for (GeneId gene : cluster) {
			unclusteredGenes.erase(gene);
		}
	}

	// Write output file
	FILE *fp = fopen(d->arguments.outputFilename.c_str(), "w");
	fprintf(fp, "Gene,Cluster\n");
	for (int i = 0; i < (int)result.clusters.size(); i++) {
		const char groupLetter = 'A' + i;
		for (const GeneId gene : result.clusters[i]) {
			const std::string geneName = d->geneRegistry.name(gene);
			fprintf(fp, "%s,%c\n", geneName.c_str(), groupLetter);
		}
	}
	// And the unclustered ones, with <null> cluster name
	for (const GeneId gene : unclusteredGenes) {
		const std::string geneName = d->geneRegistry.name(gene);
		fprintf(fp, "%s,\n", geneName.c_str());
	}
	fclose(fp);

	// Output statistics, if requested
	if (!d->arguments.statsOutputFilename.empty()) {
		FILE *fp = fopen(d->arguments.statsOutputFilename.c_str(), "w");
		fprintf(fp, "Count,Metric,Shuffled,pValue,pAdj\n");
		for (const GenomeCutter::SampleStats &stats : result.sampleStatistics) {
			fprintf(fp, "%d,%f,%f,%f,%f\n", stats.geneCount, stats.metric, stats.shuffled,stats.pValue,
					stats.adjustedPValue);
		}
		fclose(fp);
	}

	// Report time
	std::chrono::steady_clock::time_point endTime =
		std::chrono::steady_clock::now();
	int seconds =
		std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime)
			.count();
	const int minutes = seconds / 60;
	seconds %= 60;
	printf("Elapsed time: %dm %ds\n", minutes, seconds);
}

void Application::Private::initializeObjects() {
	// Load signal
	GeneSignalFactory signalFactory;
	signal = signalFactory.create(arguments.signalOptions, geneRegistry);

	// Load 3D model
	model = GenomeModel(arguments.modelOptions, geneRegistry);

	// Handshake model & signal by removing from the model nonexistent genes in signal
	const GeneSet genesInSignal = signal->geneSet();
	model.removeGenesNotExistingIn(genesInSignal, geneRegistry);

}
