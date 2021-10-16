#ifndef __APPLICATION_H__
#define __APPLICATION_H__

#include <utils/AtomBox.h>

#include "GeneSignal.h"
#include "GenomeCutter.h"
#include "GenomeModel.h"

#include <vector>

class Application {
  public:
	struct Arguments {
		GenomeModel::Options modelOptions;
		GeneSignal::Options signalOptions;
		GenomeCutter::Options cutterOptions;

		std::string outputFilename;

	};

	Application(const Arguments &arguments);
	~Application();

	void execute();

  private:
	struct Private;
	Private *d = nullptr;
};

#endif // __APPLICATION_H__