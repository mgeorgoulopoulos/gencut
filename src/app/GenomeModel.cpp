#include "GenomeModel.h"

#include "StringConstants.h"
#include <utils/Exception.h>

#include <rapidcsv.h>

void GenomeModel::Options::printHelp() const {
	printf("\t--model filename\t\tLoad 3D positions per gene from CSV "
		   "file 'filename'.\n");
}

GenomeModel::GenomeModel() {}

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