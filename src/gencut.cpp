#include <app/Application.h>

#include <app/StringConstants.h>
#include <utils/Exception.h>

#include <stdio.h>
#include <vector>

Application::Arguments parseCommandLine(const std::vector<std::string> &argv);

int main(int argc, char *argv[]) {
	try {
		// Parse command line
		std::vector<std::string> tokens;
		for (int i = 1; i < argc; i++)
			tokens.push_back(argv[i]);

		Application::Arguments arguments = parseCommandLine(tokens);

		Application app(arguments);

		app.execute();
	} catch (Exception e) {
		e.print();
		return -1;
	}

	printf("Success\n");

	return 0;
}

Application::Arguments parseCommandLine(const std::vector<std::string> &argv) {
	Application::Arguments result;

	// Parse switches
	std::vector<std::string> binaryArguments;
	for (const std::string &token : argv) {
		if (token == "--help") {
			result.displayHelp = true;
			continue;
		}
		binaryArguments.push_back(token);
	}

	if (binaryArguments.size() % 2 != 0) {
		throw(Exception({"All arguments must come in pairs like: --arg value. "
						 "You passed an odd number of command line tokens."}));
	}

	// Get pairs of --arg value
	struct NamedArgument {
		std::string name;
		std::string value;
	};

	std::vector<NamedArgument> namedArguments;

	for (int i = 0; i < binaryArguments.size(); i += 2) {
		NamedArgument arg;
		arg.name = binaryArguments[i];
		arg.value = binaryArguments[i + 1];

		if (arg.name.size() < 2 || arg.name[0] != '-' || arg.name[1] != '-') {
			throw(Exception({"Misplaced token in argument list: ", arg.name,
							 "Expected --argName"}));
		}

		// Ok. It is well formed
		arg.name = arg.name.substr(2);
		namedArguments.push_back(arg);
	}

	// Required arguments
	bool hasModel = false;
	bool hasSignal = false;
	bool hasSignalType = false;
	bool hasMetric = false;
	bool hasOutput = false;

	// Now, convert args to application format
	for (const NamedArgument &arg : namedArguments) {
		if (arg.name == kmodel) {
			result.modelOptions.filename = arg.value;
			hasModel = true;
			continue;
		}

		if (arg.name == ksignal) {
			result.signalOptions.dataFilename = arg.value;
			hasSignal = true;
			continue;
		}

		if (arg.name == ksignal_type) {
			result.signalOptions.signalType = arg.value;
			hasSignalType = true;
			continue;
		}

		if (arg.name == kmetric) {
			result.signalOptions.metric = arg.value;
			hasMetric = true;
			continue;
		}

		if (arg.name == koutput) {
			result.outputFilename = arg.value;
			hasOutput = true;
			continue;
		}

		if (arg.name == kmin_genes) {
			try {
				result.cutterOptions.minimumGeneCount = stoi(arg.value);
			} catch (...) {
				throw(Exception(
					{"Invalid value for named argument", arg.name, arg.value}));
			}
			continue;
		}

		if (arg.name == kradius) {
			try {
				result.cutterOptions.sphereRadius = stof(arg.value);
			} catch (...) {
				throw(Exception(
					{"Invalid value for named argument", arg.name, arg.value}));
			}
			continue;
		}

		if (arg.name == ksamples) {
			try {
				result.cutterOptions.randomSphereCount = stoi(arg.value);
			} catch (...) {
				throw(Exception(
					{"Invalid value for named argument", arg.name, arg.value}));
			}
			continue;
		}

		if (arg.name == kp_samples) {
			try {
				result.cutterOptions.distributionExtractionIterations =
					stoi(arg.value);
			} catch (...) {
				throw(Exception(
					{"Invalid value for named argument", arg.name, arg.value}));
			}
			continue;
		}

		if (arg.name == kpadj_threshold) {
			try {
				result.cutterOptions.pAdjThreshold = stof(arg.value);
			} catch (...) {
				throw(Exception(
					{"Invalid value for named argument", arg.name, arg.value}));
			}
			continue;
		}

		if (arg.name == koverlap_threshold) {
			try {
				result.cutterOptions.overlapThreshold = stof(arg.value);
			} catch (...) {
				throw(Exception(
					{"Invalid value for named argument", arg.name, arg.value}));
			}
			continue;
		}

		if (arg.name == ktail) {
			if (arg.value == khigh) {
				result.cutterOptions.tailSelection =
					GenomeCutter::Options::TailSelection::High;
			} else if (arg.value == klow) {
				result.cutterOptions.tailSelection =
					GenomeCutter::Options::TailSelection::Low;
			} else if (arg.value == kboth) {
				result.cutterOptions.tailSelection =
					GenomeCutter::Options::TailSelection::Both;
			} else {
				throw(Exception(
					{"Unrecognized argument for tail selection: ", arg.value}));
			}
			continue;
		}

		throw(Exception({"Unrecognized argument: ", arg.name}));
	} // end for (all argument pairs)

	if (!hasModel) {
		throw(Exception("No model file specified."));
	}

	if (!hasSignal) {
		throw(Exception("No signal file specified."));
	}

	if (!hasSignalType) {
		throw(Exception("No signal type specified."));
	}

	if (!hasMetric) {
		throw(Exception("No metric specified."));
	}

	if (!hasOutput) {
		throw(Exception("No output file specified."));
	}

	return result;
}