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
		bool displayHelp = false;

		GenomeModel::Options modelOptions;
		GeneSignal::Options signalOptions;
		GenomeCutter::Options cutterOptions;

		std::string outputFilename;
		std::string statsOutputFilename;

		void printHelp() const;
		void print() const;

	};

	Application(const Arguments &arguments);
	~Application();

	void execute();

  private:
	struct Private;
	Private *d = nullptr;
};

#endif // __APPLICATION_H__