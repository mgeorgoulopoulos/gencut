#include "GeneSignalMatrix.h"

#include "StringConstants.h"

#include <utils/Exception.h>

#include <rapidcsv.h>

#include <algorithm>
#include <functional>

struct GeneSignalMatrix::Private {
	// The huge matrix
	std::vector<double> matrix;

	int stride = 0;

	double &lookup(GeneId row, GeneId col) {
		return matrix[col + row * stride];
	}

	const double &lookup(GeneId row, GeneId col) const {
		return matrix[col + row * stride];
	}

	std::function<double(const GeneList &)> metricFunc = [](const GeneList &) {
		return 0.0;
	};

	GeneSet genes;

	// Metrics
	double average(const GeneList &genes) const;
};

#ifndef GENCUT_R_PACKAGE
GeneSignalMatrix::GeneSignalMatrix(const Options &options,
								   AtomBox &geneRegistry)
	: GeneSignal(options, geneRegistry), d(new Private()) {
	// Load CSV file
	printf("Loading matrix CSV file: %s\n", options.dataFilename.c_str());
	rapidcsv::Document doc(options.dataFilename, rapidcsv::LabelParams(0, 0));

	// Verify that we have a square matrix - ie row labels is the same set
	// as column labels
	std::set<std::string> rows;
	for (const std::string &row : doc.GetRowNames()) {
		rows.insert(row);
	}

	if (rows.empty()) {
		throw(Exception({"CSV file has zero rows:", options.dataFilename}));
	}

	std::set<std::string> columns;
	for (const std::string &col : doc.GetColumnNames()) {
		columns.insert(col);
	}

	if (rows != columns) {
		throw(Exception({"Matrix file:", options.dataFilename,
						 "Set of row names is different than set of column "
						 "names. A square matrix is required."}));
	}

	// Get set of gene Ids
	for (const std::string &gene : rows) {
		GeneId id = geneRegistry.atom(gene);
		d->genes.insert(id);
	}

	// Get greatest gene Id to calculate space required. We will flatten
	// everything to a single array.
	d->stride = *std::max_element(d->genes.begin(), d->genes.end()) + 1;

	// Calculate space required. It would be a good idea to let the user know.
	const int mb = sizeof(double) * d->stride * d->stride / 1024 / 1024;
	printf("Allocating ~%d MB of RAM for the signal matrix.\n", mb);

	// Allocate the monster and clear
	d->matrix.resize(d->stride * d->stride);
	for (int i = 0; i < d->stride * d->stride; i++)
		d->matrix[i] = 0.0;

	// Create indexes for rows and columns ids in the CSV. Loading by name is
	// far too slow.
	std::vector<std::string> csvRowNames = doc.GetRowNames();
	std::map<GeneId, int> geneIdToRowIndex;
	for (int i = 0; i < (int)csvRowNames.size(); i++) {
		const GeneId geneId = geneRegistry.atom(csvRowNames[i]);
		geneIdToRowIndex[geneId] = i;
	}
	std::vector<std::string> csvColNames = doc.GetColumnNames();
	std::map<GeneId, int> geneIdToColIndex;
	for (int i = 0; i < (int)csvColNames.size(); i++) {
		const GeneId geneId = geneRegistry.atom(csvColNames[i]);
		geneIdToColIndex[geneId] = i;
	}

	// Load cells
	for (auto rowIt = geneIdToRowIndex.begin(); rowIt != geneIdToRowIndex.end();
		 rowIt++) {
		for (auto colIt = geneIdToColIndex.begin();
			 colIt != geneIdToColIndex.end(); colIt++) {
			d->lookup(rowIt->first, colIt->first) =
				doc.GetCell<double>(rowIt->second, colIt->second);
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
#endif // ndef GENCUT_R_PACKAGE

#ifndef SKIP_R_API
GeneSignalMatrix::GeneSignalMatrix(Rcpp::NumericMatrix m, const std::string &metric, AtomBox &geneRegistry)
	: GeneSignal(GeneSignal::Options(), geneRegistry), d(new Private()) {
		
	options.metric = metric;
		
using namespace Rcpp;
		
	// Load CSV file
	printf("Loading matrix from R: rows: %d, cols: %d\n", m.nrow(), m.ncol());

	// Verify that we have a square matrix - ie row labels is the same set
	// as column labels
	std::set<std::string> rows;
	CharacterVector rnames = rownames(m);
	for (int i=0; i<rnames.size(); i++) {
		const std::string row = as<std::string>(rnames[i]);
		rows.insert(row);
	}

	if (rows.empty()) {
		throw(Exception("R matrix has zero rows"));
	}

	std::set<std::string> columns;
	CharacterVector cnames = colnames(m);
	for (int i=0; i<cnames.size(); i++) {
		const std::string col = as<std::string>(cnames[i]);
		columns.insert(col);
	}

	if (rows != columns) {
		throw(Exception("Matrix from R: Set of row names is different than set of column "
						 "names. A square matrix is required."));
	}

	// Get set of gene Ids
	for (const std::string &gene : rows) {
		GeneId id = geneRegistry.atom(gene);
		d->genes.insert(id);
	}

	// Get greatest gene Id to calculate space required. We will flatten
	// everything to a single array.
	d->stride = *std::max_element(d->genes.begin(), d->genes.end()) + 1;

	// Calculate space required. It would be a good idea to let the user know.
	const int mb = sizeof(double) * d->stride * d->stride / 1024 / 1024;
	printf("Allocating ~%d MB of RAM for the signal matrix.\n", mb);

	// Allocate the monster and clear it to zero
	d->matrix.resize(d->stride * d->stride);
	for (int i = 0; i < d->stride * d->stride; i++)
		d->matrix[i] = 0.0;

	// Create indexes for rows and columns ids in the R matrix. Loading by name is
	// far too slow.
	std::map<GeneId, int> geneIdToRowIndex;
	for (int i = 0; i < rnames.size(); i++) {
		const std::string row = as<std::string>(rnames[i]);
		const GeneId geneId = geneRegistry.atom(row);
		geneIdToRowIndex[geneId] = i;
	}
	std::map<GeneId, int> geneIdToColIndex;
	for (int i = 0; i < cnames.size(); i++) {
		const std::string col = as<std::string>(cnames[i]);
		const GeneId geneId = geneRegistry.atom(col);
		geneIdToColIndex[geneId] = i;
	}

	// Load cells
	for (auto rowIt = geneIdToRowIndex.begin(); rowIt != geneIdToRowIndex.end();
		 rowIt++) {
		for (auto colIt = geneIdToColIndex.begin();
			 colIt != geneIdToColIndex.end(); colIt++) {
			d->lookup(rowIt->first, colIt->first) =
				m(rowIt->second, colIt->second);
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
#endif // ndef SKIP_R_API

GeneSignalMatrix::~GeneSignalMatrix() { delete d; }

GeneSet GeneSignalMatrix::geneSet() const { return d->genes; }

double GeneSignalMatrix::metric(const GeneList &genes) const {
	return d->metricFunc(genes);
}

double GeneSignalMatrix::cell(GeneId a, GeneId b) const {
	return d->lookup(a, b);
}

double GeneSignalMatrix::Private::average(const GeneList &genes) const {
	if (genes.size() < 2) {
		throw(Exception("Empty or size 1 gene list provided for averaging."));
	}
	double average = 0.0;
	int pairCount = 0;
	for (int j = 0; j < (int)genes.size(); j++) {
		for (int i = 0; i < (int)genes.size(); i++) {
			if (i == j)
				continue;
			pairCount++;
			average += lookup(genes[j], genes[i]);
		}
	}
	average /= (double)pairCount;
	return average;
}