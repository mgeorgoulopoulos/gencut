#include "GenomeModel.h"

#include "StringConstants.h"
#include <utils/Exception.h>

#include <rapidcsv.h>

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
		try {
			gene.locus.chromosomeId = doc.GetCell<int>(kChromosome, i);
			gene.locus.startPosition = doc.GetCell<int>(kStart, i);
			gene.locus.endPosition = doc.GetCell<int>(kEnd, i);
		} catch (...) {
			// Do not pay much attention if these are not populated
		}

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
