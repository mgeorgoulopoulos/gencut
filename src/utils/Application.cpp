#include "Application.h"

#include "GenomeModel.h"

void Application::execute() {

	GenomeModel model;
	model.load(arguments.genomeModelFilename, geneRegistry);
}