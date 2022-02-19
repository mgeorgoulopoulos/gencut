#include "GeneSignalList.h"

#include "StringConstants.h"

#include <utils/Exception.h>

#include <rapidcsv.h>

#include <algorithm>
#include <functional>

struct GeneSignalList::Private {
	AtomBox &geneRegistry;
	std::map<GeneId, double> geneList;
	GeneSet genes;

	Private(AtomBox &geneRegistry) : geneRegistry(geneRegistry) {}

	std::function<double(const GeneList &)> metricFunc = [](const GeneList &) {
		return 0.0;
	};

	// Returns true for metrics that require factor input (ex: 'entropy').
	static bool isFactorMetric(const std::string &metric) {
		if (metric == kaverage)
			return false;

		return false;
	}

	// Metrics
	double average(const GeneList &genes) const;
};

GeneSignalList::GeneSignalList(const Options &options, AtomBox &geneRegistry)
	: GeneSignal(options, geneRegistry), d(new Private(geneRegistry)) {
	// Load CSV file
	printf("Loading gene list CSV file: %s\n", options.dataFilename.c_str());
	rapidcsv::Document doc(
		options.dataFilename,
		rapidcsv::LabelParams(0 /* First row = column naems */));

	// Verify that the CSV has the correct columns
	const std::vector<std::string> columnNames = doc.GetColumnNames();
	const std::vector<std::string> expectedColumnNames = {kgene, kvalue};

	if (columnNames.size() != expectedColumnNames.size()) {
		throw(Exception({options.dataFilename, "2 columns expected"}));
	}

	for (int i = 0; i < (int)columnNames.size(); i++) {
		if (columnNames[i] != expectedColumnNames[i]) {
			throw(Exception({"Settings file:", options.dataFilename,
							 "Expected column:", expectedColumnNames[i],
							 "Actual column:", columnNames[i]}));
		}
	}

	// Load data
	const bool factorData = Private::isFactorMetric(options.metric);
	std::map<std::string, int> factors;
	for (int rowId = 0; rowId < (int)doc.GetRowCount(); rowId++) {
		const std::string gene = doc.GetCell<std::string>(0, rowId);
		const GeneId id = geneRegistry.atom(gene);
		d->genes.insert(id);

		if (factorData) {
			int factored = 0;
			const std::string data = doc.GetCell<std::string>(1, rowId);
			auto it = factors.find(data);
			if (it == factors.end()) {
				factored = (int)factors.size();
				factors[data] = factored;
			} else {
				factored = it->second;
			}
			d->geneList[id] = (double)factored;
		} else {
			d->geneList[id] = doc.GetCell<double>(1, rowId);
		}
	}

	// Finally, setup metric
	if (options.metric == kaverage) {
		d->metricFunc = [&](const GeneList &genes) {
			return d->average(genes);
		};
	} else {
		throw(Exception(
			{"Matrix signal does not support this metric:", options.metric}));
	}
}

GeneSignalList::~GeneSignalList() { delete d; }

GeneSet GeneSignalList::geneSet() const { return d->genes; }

double GeneSignalList::metric(const GeneList &genes) const {
	return d->metricFunc(genes);
}

double GeneSignalList::Private::average(const GeneList &genes) const {
	if (genes.empty()) {
		throw(Exception("Empty gene list provided for averaging."));
	}
	double average = 0.0;
	for (const GeneId &geneId : genes) {
		const auto it = geneList.find(geneId);
		if (it == geneList.end()) {
			throw(Exception(
				{"average: gene is not contained in the provided list: ",
				 geneRegistry.name(geneId)}));
		}
		average += it->second;
	}
	average /= (double)genes.size();
	return average;
}