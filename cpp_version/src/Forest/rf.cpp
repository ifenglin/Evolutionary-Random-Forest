/*-------------------------------------------------------------------------------
 This file is part of Ranger.

 Ranger is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Ranger is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with Ranger. If not, see <http://www.gnu.org/licenses/>.

 Written by:

 Marvin N. Wright
 Institut für Medizinische Biometrie und Statistik
 Universität zu Lübeck
 Ratzeburger Allee 160
 23562 Lübeck
 Germany

 http://www.imbs-luebeck.de
 #-------------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <iostream>
#include "rf.h"

using namespace std;


Forest* rf(int argc, char **argv, std::vector<genotype> &genes, bool reload_data) {
	//cout << "Start growing forest: " << endl;
	//for (int i = 0; i < argc ; i++) {
	//	cout << argv[i] << " ";
	//}
	//cout << endl;
  ArgumentHandler arg_handler(argc, argv);
  Forest* forest = 0;
  static Data* data = 0;
  static std::vector<double>* case_weights = 0;

  try {

    // Handle command line arguments
    if (arg_handler.processArguments() != 0) {
		cout << "Cannot handle command line arguments." << endl;
      return 0;
    }
	//cout << "processed" << endl << "About to check arguments...";
    arg_handler.checkArguments();
	//cout << "arguments checked" << endl;
    // Create forest object
    /*switch (arg_handler.treetype) {
    case TREE_CLASSIFICATION:
      if (arg_handler.probability) {
        forest = new ForestProbability;
      } else {
        forest = new ForestClassification;
      }
      break;
    case TREE_REGRESSION:
      forest = new ForestRegression;
      break;
    case TREE_SURVIVAL:
      forest = new ForestSurvival;
      break;
    case TREE_PROBABILITY:
      forest = new ForestProbability;
      break;
    }*/
	forest = new ForestClassification;

    // Verbose output to logfile if non-verbose mode
    static std::ostream* verbose_out;
    if (arg_handler.verbose) {
      verbose_out = &std::cout;
    } else {
      static std::ofstream* logfile = new std::ofstream();
	  if (verbose_out != logfile) {
		  logfile->open(arg_handler.outprefix + "_log.txt");
		  if (!logfile->good()) {
			  throw std::runtime_error("Could not write to logfile.");
		  }
		  verbose_out = logfile;
	  }
    }

    // Call Ranger
    *verbose_out << "Starting Ranger." << std::endl;
	if (reload_data) {
		data = forest->loadData(arg_handler.file);
		*verbose_out << "Load data...";
	}
	*verbose_out << "Set data...";
	forest->setData(data);
	*verbose_out << "done." << endl;
	if (reload_data) {
		delete case_weights;
		*verbose_out << "Load case weights file...";
		case_weights = forest->loadCaseWeights(arg_handler.caseweights);
	}
	*verbose_out << "Set case weights...";
	forest->setCaseWeights(case_weights);
	*verbose_out << "done." << std::endl << "Set genes...";
	forest->setGenes(genes);
	*verbose_out << "done." << std::endl << "Initialize forest...";
    forest->initCpp(arg_handler.depvarname, arg_handler.memmode, arg_handler.file, arg_handler.mtry,
        arg_handler.outprefix, arg_handler.ntree, verbose_out, arg_handler.seed, arg_handler.nthreads,
        arg_handler.predict, arg_handler.impmeasure, arg_handler.targetpartitionsize, arg_handler.splitweights,
        arg_handler.alwayssplitvars, arg_handler.statusvarname, arg_handler.replace, arg_handler.catvars,
        arg_handler.savemem, arg_handler.splitrule, arg_handler.caseweights, arg_handler.predall, arg_handler.fraction,
        arg_handler.alpha, arg_handler.minprop, arg_handler.holdout, arg_handler.predictiontype,
        arg_handler.randomsplits);
	*verbose_out << "forest initilized." << std::endl << "Run forest...";
    forest->run(true);
	*verbose_out << "finished." << std::endl;
    if (arg_handler.write) {
      forest->saveToFile();
    }
    forest->writeOutput();
    *verbose_out << "Finished Ranger." << std::endl;

    //delete forest;
  } catch (std::exception& e) {
    std::cerr << "Error: " << e.what() << " Ranger will EXIT now." << std::endl;
    delete forest;
	forest = NULL;
  }

  return forest;
}
