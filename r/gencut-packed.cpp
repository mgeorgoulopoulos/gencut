// ---------------------------------
// Start file: includes.cpp
// ---------------------------------

#include <omp.h>

#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <cmath>

#define GENCUT_R_PACKAGE

// [[Rcpp::plugins(openmp)]]

#ifndef SKIP_R_API
#include <Rcpp.h>
#define printf Rprintf
#endif
// ---------------------------------
// Start file: ../src/utils/AtomBox.h
// ---------------------------------

#ifndef __ATOM_BOX_H__
#define __ATOM_BOX_H__


using Atom = int;

class AtomBox : public std::map<std::string, Atom> {
  public:
	bool contains(const std::string &geneName) const {
		return find(geneName) != end();
	}

	Atom atom(const std::string &geneName) {
		auto it = find(geneName);
		if (it != end())
			return it->second;

		// Not found - insert
		const Atom atom = (int)size() + 1;
		(*this)[geneName] = atom;

		return atom;
	}

	// Too bored to do this properly now
	std::string name(Atom atom) {
		for (auto it = begin(); it != end(); it++) {
			if (it->second == atom)
				return it->first;
		}
		return "";
	}
};

#endif // __ATOM_BOX_H__
// ---------------------------------
// Start file: ../src/utils/Exception.h
// ---------------------------------

#ifndef __EXCEPTION_H__
#define __EXCEPTION_H__


class Exception {
public:
	using StringList = std::vector<std::string>;
	
	Exception() {}

	Exception(const std::string &message) : Exception(StringList{message}) {}
	
	Exception(const StringList &strings) : strings(strings) {}
	
	void print() {
		printf("Exception details:\n");
		for (const std::string &s : strings) {
			printf("\t%s\n", s.c_str());			
		}
	}

private:
	StringList strings;
	
};

#endif // __EXCEPTION_H__
// ---------------------------------
// Start file: ../src/utils/Vec3D.h
// ---------------------------------

/*
Copyright 2021 Michael Georgoulopoulos

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/*
Vector class that represents points in 3D space. Initially I used GLM, but I
decided it was not worth the dependency.
*/

#ifndef _VEC_3D_H_
#define _VEC_3D_H_


class Vec3D {
  public:
	Vec3D() {}
	Vec3D(double x, double y, double z) : x(x), y(y), z(z) {}

	double x = 0.0;
	double y = 0.0;
	double z = 0.0;

	// Commonly used math operators
	Vec3D operator+(const Vec3D &other) const {
		return Vec3D(x + other.x, y + other.y, z + other.z);
	}

	void operator+=(const Vec3D &other) { *this = *this + other; }

	Vec3D operator-(const Vec3D &other) const {
		return Vec3D(x - other.x, y - other.y, z - other.z);
	}

	void operator-=(const Vec3D &other) { *this = *this - other; }

	Vec3D operator*(double multiplier) const {
		return Vec3D(x * multiplier, y * multiplier, z * multiplier);
	}

	void operator*=(double multiplier) { *this = *this * multiplier; }

	Vec3D operator/(double divisor) const {
		const double multiplier = 1.0 / divisor;
		return *this * multiplier;
	}

	void operator/=(double divisor) { *this = *this / divisor; }

	// Length (magnitude) of the vector
	double length() const { return sqrt(x * x + y * y + z * z); }

	// Euclidean distance
	static double distance(const Vec3D &a, const Vec3D &b) {
		const Vec3D displacement = b - a;
		return displacement.length();
	}

	// Linear interpolation of scalars
	static double mix(double a, double b, double t) {
		const double omt = 1.0 - t;
		return a * omt + b * t;
	}

	static Vec3D mix(const Vec3D &a, const Vec3D &b, double t) {
		return Vec3D(mix(a.x, b.x, t), mix(a.y, b.y, t), mix(a.z, b.z, t));
	}
};

#endif // _VEC_3D_H_

// ---------------------------------
// Start file: ../src/app/StringConstants.cpp
// ---------------------------------


const char *kaverage = "average";
const char *kboth = "both";
const char *kChromosome = "Chromosome";
const char *kEnd = "End";
const char *kentropy = "entropy";
const char *kGene = "Gene";
const char *khigh = "high";
const char *kkey = "key";
const char *klist = "list";
const char *klow = "low";
const char *kmatrix = "matrix";
const char *kmetric = "metric";
const char *kmin_genes = "min-genes";
const char *kmodel = "model";
const char *kp_samples = "p-samples";
const char *ksamples = "samples";
const char *ksignal = "signal";
const char *ksignal_type = "signal-type";
const char *kStart = "Start";
const char *kstats_output = "stats-output";
const char *ktail = "tail";
const char *koutput = "output";
const char *koverlap_threshold = "overlap-threshold";
const char *kpadj_threshold = "padj-threshold";
const char *kradius = "radius";
const char *ksettings = "settings";
const char *kValue = "Value";
const char *kx = "x";
const char *ky = "y";
const char *kz = "z";


// ---------------------------------
// Start file: ../src/app/GeneCollections.h
// ---------------------------------

#ifndef __GENE_COLLECTIONS_H__
#define __GENE_COLLECTIONS_H__



using GeneId = Atom;

class GeneSet : public std::set<GeneId> {
  public:
	bool contains(GeneId g) { return find(g) != end();
	}
};

class GeneList : public std::vector<GeneId> {};

#endif // __GENE_COLLECTIONS_H__
// ---------------------------------
// Start file: ../src/app/GeneSignal.h
// ---------------------------------

#ifndef __GENE_SIGNAL_H__
#define __GENE_SIGNAL_H__



// Abstract class. Calculates metric for sets of genes.
class GeneSignal {
  public:
	struct Options {
		std::string dataFilename;
		std::string signalType;
		std::string metric;

		void print() const;
		void printHelp() const;
	};
	GeneSignal(const Options &options, AtomBox &geneRegistry);
	virtual ~GeneSignal();

	virtual GeneSet geneSet() const = 0;

	virtual double metric(const GeneList &genes) const = 0;

  protected:
	Options options;
	AtomBox &geneRegistry;
};

using GeneSignalPtr = std::shared_ptr<GeneSignal>;

#endif // __GENE_SIGNAL_H__
// ---------------------------------
// Start file: ../src/app/GeneSignal.cpp
// ---------------------------------


void GeneSignal::Options::print() const {
	printf("\tSignal type: %s\n", signalType.c_str());
	printf("\tMetric: %s\n", metric.c_str());
}

void GeneSignal::Options::printHelp() const {
	printf("\t--signal filename\t\tLoad gene signal from CSV file 'filename'.\n");
	printf("\t--signal-type type\t\tSpecify type of signal. Currently only 'matrix' is valid.\n");
	printf("\t--metric function\t\tSpecify function to be applied to sets of genes, ex: 'average'.\n");
}

GeneSignal::GeneSignal(const Options &options, AtomBox &geneRegistry)
	: options(options), geneRegistry(geneRegistry) {}
GeneSignal::~GeneSignal() {}
// ---------------------------------
// Start file: ../src/app/GeneSignalMatrix.h
// ---------------------------------

#ifndef __GENE_SIGNAL_MATRIX_H__
#define __GENE_SIGNAL_MATRIX_H__


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
// ---------------------------------
// Start file: ../src/app/GeneSignalMatrix.cpp
// ---------------------------------






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
// ---------------------------------
// Start file: ../src/app/GenomeModel.h
// ---------------------------------

#ifndef __GENOME_MODEL_H__
#define __GENOME_MODEL_H__




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
// ---------------------------------
// Start file: ../src/app/GenomeModel.cpp
// ---------------------------------




void GenomeModel::Options::printHelp() const {
	printf("\t--model filename\t\tLoad 3D positions per gene from CSV "
		   "file 'filename'.\n");
}

GenomeModel::GenomeModel() {}

#ifndef GENCUT_R_PACKAGE

GenomeModel::GenomeModel(const Options &options, AtomBox &geneRegistry) {
	printf("Loading genome 3D model from CSV file: %s\n",
		   options.filename.c_str());
	rapidcsv::Document doc(options.filename);

	// Verify that the CSV has the correct columns
	const std::vector<std::string> columnNames = doc.GetColumnNames();
	const std::vector<std::string> expectedColumnNames = {
		kGene, kChromosome, kStart, kEnd, kx, ky, kz};

	if (columnNames.size() != expectedColumnNames.size()) {
		throw(Exception({options.filename, "7 columns expected"}));
	}

	for (int i = 0; i < (int)columnNames.size(); i++) {
		if (columnNames[i] != expectedColumnNames[i]) {
			throw(Exception({"Genome model file:", options.filename,
							 "Expected column:", expectedColumnNames[i],
							 "Actual column:", columnNames[i]}));
		}
	}

	// Load records
	genes.resize(0);
	genes.reserve(doc.GetRowCount());
	for (int i = 0; i < (int)doc.GetRowCount(); i++) {
		const std::string &geneName = doc.GetCell<std::string>(kGene, i);

		Gene gene;
		gene.geneId = geneRegistry.atom(geneName);
		gene.locus.chromosomeId = doc.GetCell<int>(kChromosome, i);
		gene.locus.startPosition = doc.GetCell<int>(kStart, i);
		gene.locus.endPosition = doc.GetCell<int>(kEnd, i);
		gene.position.x = doc.GetCell<double>(kx, i);
		gene.position.y = doc.GetCell<double>(ky, i);
		gene.position.z = doc.GetCell<double>(kz, i);

		genes.push_back(gene);
	}

	printf("%d genes loaded\n", (int)genes.size());
}

#endif // ndef GENCUT_R_PACKAGE

#ifndef SKIP_R_API
GenomeModel::GenomeModel(Rcpp::NumericMatrix m, AtomBox &geneRegistry) {

	using namespace Rcpp;
	printf("Loading genome 3D model from R matrix: rows: %d, cols: %d\n",
		   m.nrow(), m.ncol());

	// Verify that the matrix has the correct columns

	if (m.ncol() != 3) {
		throw(Exception("R matrix does not have 3 columns. For a genome model "
						"3 columns (x,y,z) are required."));
	}

	// Load records
	CharacterVector rnames = rownames(m);
	genes.resize(0);
	genes.reserve(m.nrow());
	for (int i = 0; i < m.nrow(); i++) {
		const std::string &geneName = as<std::string>(rnames[i]);

		Gene gene;
		gene.geneId = geneRegistry.atom(geneName);

		// Leave chromosome / start / end empty for now
		// gene.locus.chromosomeId = doc.GetCell<int>(kChromosome, i);
		// gene.locus.startPosition = doc.GetCell<int>(kStart, i);
		// gene.locus.endPosition = doc.GetCell<int>(kEnd, i);
		gene.position.x = m(i, 0);
		gene.position.y = m(i, 1);
		gene.position.z = m(i, 2);

		genes.push_back(gene);
	}

	printf("%d genes loaded\n", (int)genes.size());
}
#endif // ndef SKIP_R_API

GeneSet GenomeModel::geneSet() const {
	GeneSet result;
	for (const Gene &gene : genes) {
		result.insert(gene.geneId);
	}
	return result;
}

GeneList GenomeModel::geneList() const {
	GeneList result;
	for (const Gene &gene : genes) {
		result.push_back(gene.geneId);
	}
	return result;
}

void GenomeModel::removeGenes(const GeneSet &genesToRemove) {
	std::vector<Gene> newGenes;
	for (const Gene &gene : genes) {
		if (genesToRemove.find(gene.geneId) != genesToRemove.end()) {
			continue;
		}
		newGenes.push_back(gene);
	}
	genes = newGenes;
}

GeneList GenomeModel::sphereSample(const Vec3D &center, double radius) const {
	GeneList result;
	for (const Gene &gene : genes) {
		const double distance = Vec3D::distance(center, gene.position);
		if (distance <= radius) {
			result.push_back(gene.geneId);
		}
	}
	return result;
}

Vec3D GenomeModel::minBound() const {
	if (genes.empty()) {
		throw(Exception("Genome model is empty."));
	}
	Vec3D result = genes.front().position;
	for (const Gene &gene : genes) {
		result.x = std::min(result.x, gene.position.x);
		result.y = std::min(result.y, gene.position.y);
		result.z = std::min(result.z, gene.position.z);
	}
	return result;
}

Vec3D GenomeModel::maxBound() const {
	if (genes.empty()) {
		throw(Exception("Genome model is empty."));
	}
	Vec3D result = genes.front().position;
	for (const Gene &gene : genes) {
		result.x = std::max(result.x, gene.position.x);
		result.y = std::max(result.y, gene.position.y);
		result.z = std::max(result.z, gene.position.z);
	}
	return result;
}

void GenomeModel::removeGenesNotExistingIn(const GeneSet &genes, AtomBox &geneRegistry) {
	// Negotiate gene lists between signal and model. We will remove nonexistent
	// genes in signal from model.
	GeneSet genesInModel = geneSet();
	GeneSet difference;
	std::set_difference(genesInModel.begin(), genesInModel.end(),
						genes.begin(), genes.end(),
						std::inserter(difference, difference.begin()));

	if (!difference.empty()) {
		printf("Removing %d genes from the 3D model because they don't exist "
			   "in the signal.\n",
			   difference.size());
		for (GeneId gene : difference) {
			const std::string name = geneRegistry.name(gene);
			printf("%s ", name.c_str());
		}
		printf("\n");
		removeGenes(difference);
		printf("3D model remains with %d genes\n", geneSet().size());
	}
}

// ---------------------------------
// Start file: ../src/app/GenomeCutter.h
// ---------------------------------

#ifndef __GENOME_CUTTER_H__
#define __GENOME_CUTTER_H__




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

		static TailSelection tailSelectionFromString(const std::string &s);

		void print() const;
		void printHelp() const;
	};

	GenomeCutter(const GenomeModel &model, const GeneSignalPtr signal,
				 const Options &options);
	virtual ~GenomeCutter();

	// Statistics of a sphere-sample.
	struct SampleStats {
		double metric = 0.0;
		double shuffled = 0.0; // value of metric coming from shuffled gene set of the same size
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
// ---------------------------------
// Start file: ../src/app/GenomeCutter.cpp
// ---------------------------------






GenomeCutter::Options::TailSelection
GenomeCutter::Options::tailSelectionFromString(const std::string &s) {
	if (s == khigh) {
		return GenomeCutter::Options::TailSelection::High;
	} else if (s == klow) {
		return GenomeCutter::Options::TailSelection::Low;
	} else if (s == kboth) {
		return GenomeCutter::Options::TailSelection::Both;
	}

	throw(Exception({"Unrecognized argument for tail selection: ", s}));
}

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
	double shuffled = 0.0; // a value of the metric from a shuffled gene set of the same size

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

GenomeCutter::Result GenomeCutter::cut() const {
	Result result;

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
		const std::vector<double> &randomizedSample =
			geneCountToDistribution[(int)workUnit.genes.size()];
		if (randomizedSample.empty()) {
			throw(Exception("Empty randomized sample"));
		}
		workUnit.metric = d->signal->metric(workUnit.genes);
		workUnit.pValue = d->calculatePValue(workUnit.metric, randomizedSample);
		workUnit.shuffled = randomizedSample[rand() % randomizedSample.size()];
	}
	printf("Done\n");

	// Adjust p-values
	printf("Adjusting p-values using Benjamini-Hochberg method... ");
	benjamini(workUnits);
	printf("Done.\n");

	// Export stats
	result.sampleStatistics.reserve(workUnits.size());
	for (const WorkUnit &workUnit : workUnits) {
		SampleStats s;
		s.metric = workUnit.metric;
		s.shuffled = workUnit.shuffled;
		s.pValue = workUnit.pValue;
		s.adjustedPValue = workUnit.adjustedPValue;
		result.sampleStatistics.push_back(s);
	}

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

	if (workUnits.empty()) {
		printf("No significant samples found");
	}

	// Convert work units to clusters of genes
	printf("\n");
	printf("Hierarchical clustering. Using threshold of %.02f%% overlap ratio "
		   "to consider clusters as distinct ... ",
		   d->options.overlapThreshold * 100.0);
	double maximumOverlapRatio = 0.0;
	result.clusters = clusterByGeneOverlap(
		workUnits, d->options.overlapThreshold, &maximumOverlapRatio);
	printf("Done.\n");
	printf("Stopping clustering with %d clusters, %.02f%% maximum gene "
		   "overlap.\n",
		   result.clusters.size(), maximumOverlapRatio * 100.0);

	// The returned clusters are fairly defined. Sometimes we get a different
	// number of clusters because this is a random sampling, but mostly we get
	// the same. Ordering the clusters by size.
	std::sort(
		result.clusters.begin(), result.clusters.end(),
		[](const GeneSet &a, const GeneSet &b) { return a.size() < b.size(); });

	// Cleanup duplicate genes - let the smaller cluster retain all its genes
	GeneSet genesUsed;
	for (GeneSet &cluster : result.clusters) {
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
	for (int i = 0; i < (int)result.clusters.size(); i++) {
		printf("\tCluster %c: %d genes\n", 'A' + i, result.clusters[i].size());
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
// ---------------------------------
// Start file: rapi.cpp
// ---------------------------------

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

