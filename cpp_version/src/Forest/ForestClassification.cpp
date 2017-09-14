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

#include <unordered_map>
#include <algorithm>
#include <iterator>
#include <random>
#include <stdexcept>
#include <cmath>
#include <string>

#include "utility.h"
#include "ForestClassification.h"
#include "TreeClassification.h"
#include "Data.h"

ForestClassification::ForestClassification() {
}

ForestClassification::~ForestClassification() {
}

void ForestClassification::loadForest(size_t dependent_varID, size_t num_trees,
    std::vector<std::vector<std::vector<size_t>> >& forest_child_nodeIDs,
    std::vector<std::vector<size_t>>& forest_split_varIDs, std::vector<std::vector<double>>& forest_split_values,
    std::vector<double>& class_values, std::vector<bool>& is_ordered_variable) {

  this->dependent_varID = dependent_varID;
  this->num_trees = num_trees;
  this->class_values = class_values;
  this->is_ordered_variable = is_ordered_variable;

  // Create trees
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    Tree* tree = new TreeClassification(forest_child_nodeIDs[i], forest_split_varIDs[i], forest_split_values[i],
        &this->class_values, &response_classIDs, &this->is_ordered_variable);
    trees.push_back(tree);
  }

  // Create thread ranges
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);
}

void ForestClassification::initInternal(std::string status_variable_name) {

  // If mtry not set, use floored square root of number of independent variables.
  if (mtry == 0) {
    unsigned long temp = sqrt((double) (num_variables - 1));
    mtry = std::max((unsigned long) 1, temp);
  }

  // Set minimal node size
  if (min_node_size == 0) {
    min_node_size = DEFAULT_MIN_NODE_SIZE_CLASSIFICATION;
  }

  // Create class_values and response_classIDs
  if (!prediction_mode) {
    for (size_t i = 0; i < num_samples; ++i) {
      double value = data->get(i, dependent_varID);

      // If classID is already in class_values, use ID. Else create a new one.
      uint classID = find(class_values.begin(), class_values.end(), value) - class_values.begin();
      if (classID == class_values.size()) {
        class_values.push_back(value);
      }
      response_classIDs.push_back(classID);
    }
  }

  // Sort data if memory saving mode
  if (!memory_saving_splitting) {
    data->sort();
  }

  // initialize prediction error of each tree
  prediction_error_each_tree = std::vector<double>(num_trees);
}

void ForestClassification::growInternal() {
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    trees.push_back(new TreeClassification(&class_values, &response_classIDs));
  }
}

void ForestClassification::predictInternal() {

  size_t num_prediction_samples = data->getNumRows();
  if (predict_all || prediction_type == TERMINALNODES) {
    predictions = std::vector<std::vector<std::unordered_map<size_t, unsigned>>>(1, std::vector<std::unordered_map<size_t, unsigned>>(num_prediction_samples, std::unordered_map<size_t, unsigned>(num_trees)));
  } else {
    predictions = std::vector<std::vector<std::unordered_map<size_t, unsigned>>>(1, std::vector<std::unordered_map<size_t, unsigned>>(1, std::unordered_map<size_t, unsigned>(num_prediction_samples)));
  }

  // For all samples get tree predictions
  for (size_t sample_idx = 0; sample_idx < num_prediction_samples; ++sample_idx) {

    if (predict_all || prediction_type == TERMINALNODES) {
      // Get all tree predictions
      for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
        if (prediction_type == TERMINALNODES) {
          predictions[0][sample_idx][tree_idx] = ((TreeClassification*) trees[tree_idx])->getPredictionTerminalNodeID(
              sample_idx);
        } else {
          predictions[0][sample_idx][tree_idx] = ((TreeClassification*) trees[tree_idx])->getPrediction(sample_idx);
        }
      }
    } else {
      // Count classes over trees and save class with maximum count
      std::unordered_map<double, size_t> class_count;
      for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
        double value = ((TreeClassification*) trees[tree_idx])->getPrediction(sample_idx);
        ++class_count[value];
      }
      predictions[0][0][sample_idx] = mostFrequentValue(class_count, random_number_generator);
    }

  }
}

void ForestClassification::computePredictionErrorInternal() {

	// Use a fraction of the samples
	size_t num_samples_downsized = (size_t)(num_samples/num_trees);

	std::vector<size_t> sampleIDs;
	std::uniform_int_distribution<size_t> unif_dist(0, num_samples - 1);

	sampleIDs.reserve(num_samples_downsized);
	for (size_t s = 0; s < num_samples_downsized; ++s) {
		size_t draw = unif_dist(random_number_generator);
		sampleIDs[s] = draw;
	}

  // Class counts for samples
  std::unordered_map<size_t, std::unordered_map<double, size_t>> class_counts;
  class_counts.reserve(num_samples_downsized);
  // For each tree loop over OOB samples and count classes
  // record predictions from each tree
  *verbose_out << "-Recording predictions..";
  //predictions_each_tree = std::vector<std::vector<std::vector<double>>>(1, std::vector<std::vector<double>>(num_trees, std::vector<double>(num_samples)));
  std::vector<size_t> num_missclassifications_each_tree = std::vector<size_t>(num_trees, 0);
  std::vector<size_t> num_predictions_each_tree = std::vector<size_t>(num_trees, 0);
  //prediction_error_each_tree = std::vector<double>(num_trees);
  try {
	  for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
		  for (size_t sample_idx = 0; sample_idx < num_samples_downsized; ++sample_idx) {
			  size_t sampleID = sampleIDs[sample_idx];
			  if (trees[tree_idx]->isOobSample(sampleID)) {
				  size_t sample_idx_tree = trees[tree_idx]->findOobSample(sampleID);
				  double value = ((TreeClassification*)trees[tree_idx])->getPrediction(sample_idx_tree);
				  double real_value = data->get(sample_idx, dependent_varID);
			      ++class_counts[sampleID][value];
				  //predictions_each_tree[0][tree_idx][sampleID] = value
				  ++num_predictions_each_tree[tree_idx];
				  if (value != real_value) {
					  ++num_missclassifications_each_tree[tree_idx];
				  }
			  }
		  }
		  prediction_error_each_tree[tree_idx] = (double)num_missclassifications_each_tree[tree_idx] / (double)num_predictions_each_tree[tree_idx];
	  }
  }  catch (const std::bad_alloc &e) {
	  std::cout << "Allocation failed when recording prediction of OOB sample: " << e.what() << ". Skipped" << std::endl;
  }

  // Compute majority vote for each sample
  *verbose_out << "-Computing majority..";
  predictions = std::vector<std::vector<std::unordered_map<size_t, unsigned>>>(1, std::vector<std::unordered_map<size_t, unsigned>>(1, std::unordered_map<size_t, unsigned>(num_samples_downsized)));;
  try {
	  for (size_t i = 0; i < num_samples_downsized; ++i) {
		  size_t sampleID = sampleIDs[i];
		  if (!class_counts[sampleID].empty()) {
			  predictions[0][0][sampleID] = mostFrequentValue(class_counts[sampleID], random_number_generator);
		  }
		  else {
        // for double type prediction:
			  //predictions[0][0][i] = NAN;
        // for unsigned type prediction:
			  predictions[0][0][sampleID] = 0;
		  }
	  }
  } catch (const std::bad_alloc &e) {
	  std::cout << "Allocation failed when compute majority vote for prediction: " << e.what() << ". Skipped" << std::endl;
  }

  // Compare predictions with true data
  *verbose_out << "-Comparing predictions..";
  size_t num_missclassifications = 0;
  size_t num_predictions = 0;
  try {
	  double predicted_value, real_value;
	  for (size_t i = 0; i < num_samples_downsized; ++i) {
		  size_t sampleID = sampleIDs[i];
		  predicted_value = predictions[0][0][sampleID];
      // for double type prediction:
		  //if (!std::isnan(predicted_value)) {
      // for unsigned type prediction:
		  if(predicted_value != 0) {
			  ++num_predictions;
			  real_value = data->get(sampleID, dependent_varID);
			  if (predicted_value != real_value) {
				  ++num_missclassifications;
			  }
			  ++classification_table[std::make_pair(real_value, predicted_value)];
			  // accumulate errors for each tree
			  //for ( tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
				 // //predicted_value = predictions_each_tree[0][tree_idx][i];
				 // predicted_value = ((TreeClassification*)trees[tree_idx])->getPrediction(i);
				 // if (predicted_value != 0 && predicted_value != real_value) {
					//  ++num_missclassifications_each_tree[tree_idx];
				 // }
			  //}
		  }
	  }
  } catch (const std::bad_alloc &e) {
	  std::cout << "Allocation failed when comparing prediction with sample: " << e.what() << ". Skipped" << std::endl;
  }

  overall_prediction_error = (double) num_missclassifications / (double) num_predictions;

  //calculate error for each tree
  //*verbose_out << "-Calculating error..";
  //prediction_error_each_tree = std::vector<double>(num_trees);
  //try {
	 // for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
		//  prediction_error_each_tree[tree_idx] = (double)num_missclassifications_each_tree[tree_idx] / (double)trees[tree_idx]->getNumSamplesOob();
	 // }
  //} catch (const std::bad_alloc &e) {
	 // std::cout << "Allocation failed when calculating error for tree: " << e.what() << ". Skipped" << std::endl;
  //}

  // calculate (directed) correlation
  *verbose_out << "-Calculating correlation..";
  try {
	  correlation_each_tree = std::vector<std::vector<std::vector<double>>>(1, std::vector<std::vector<double>>(num_trees, std::vector<double>(num_trees)));
	  //std::vector<size_t> intersected_oob_sampleIDs;
	  std::vector<size_t> i_oob_sampleIDs;
	  std::vector<size_t> j_oob_sampleIDs;
	  for ( size_t i = 0; i < num_trees; ++i) {
		  i_oob_sampleIDs = trees[i]->getOobSampleIDs();
		  for ( size_t j = i + 1; j < num_trees; ++j) {
			size_t correlated_count = 0;
			size_t uncorrelated_count = 0;
			j_oob_sampleIDs = trees[j]->getOobSampleIDs();
			//intersected_oob_sampleIDs.clear();
			//std::set_intersection(i_oob_sampleIDs.begin(), i_oob_sampleIDs.end(), j_oob_sampleIDs.begin(), j_oob_sampleIDs.end(), std::back_inserter(intersected_oob_sampleIDs));
			for ( size_t sample_idx = 0; sample_idx < num_samples_downsized; ++sample_idx) {
				size_t sampleID = sampleIDs[sample_idx];
				if (trees[i]->isOobSample(sampleID) && trees[j]->isOobSample(sampleID)) {
				//if (std::find(intersected_oob_sampleIDs.begin(), intersected_oob_sampleIDs.end(), sampleID) != intersected_oob_sampleIDs.end()) {
					size_t true_value = data->get(sampleID, dependent_varID);
					size_t sample_idx_tree_i = trees[i]->findOobSample(sampleID);
					size_t sample_idx_tree_j = trees[j]->findOobSample(sampleID);
					double predicted_value_i = ((TreeClassification*)trees[i])->getPrediction(sample_idx_tree_i);
					double predicted_value_j = ((TreeClassification*)trees[j])->getPrediction(sample_idx_tree_j);
					if (predicted_value_i != predicted_value_j) {
						++uncorrelated_count;
					}
					else {
						if (predicted_value_i != true_value) {
							++correlated_count;
						}
					}
				}
			}
			double correlation_rate = (double)correlated_count / (double)(correlated_count + uncorrelated_count);
			// set correlation only to the tree that has a higher prediction error
			if (prediction_error_each_tree[i] < prediction_error_each_tree[j]) {
				correlation_each_tree[0][i][j] = 0;
				correlation_each_tree[0][j][i] = correlation_rate;
			}
			else {
				correlation_each_tree[0][i][j] = correlation_rate;
				correlation_each_tree[0][j][i] = 0;
			}
		  }
	  }
  } catch (const std::bad_alloc &e) {
	  std::cout << "Allocation failed when calculating correation between tree: " << e.what() << ". Skipped" << std::endl;
  }
}

// #nocov start
void ForestClassification::writeOutputInternal() {
  *verbose_out << "Tree type:                         " << "Classification" << std::endl;
}

void ForestClassification::writeConfusionFile() {

  // Open confusion file for writing
  std::string filename = output_prefix + ".confusion";
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to confusion file: " + filename + ".");
  }

  // Write confusion to file
  outfile << "Overall OOB prediction error (Fraction missclassified): " << overall_prediction_error << std::endl;
  outfile << std::endl;
  outfile << "Class specific prediction errors:" << std::endl;
  outfile << "           ";
  for (auto& class_value : class_values) {
    outfile << "     " << class_value;
  }
  outfile << std::endl;
  for (auto& predicted_value : class_values) {
    outfile << "predicted " << predicted_value << "     ";
    for (auto& real_value : class_values) {
      size_t value = classification_table[std::make_pair(real_value, predicted_value)];
      outfile << value;
      if (value < 10) {
        outfile << "     ";
      } else if (value < 100) {
        outfile << "    ";
      } else if (value < 1000) {
        outfile << "   ";
      } else if (value < 10000) {
        outfile << "  ";
      } else if (value < 100000) {
        outfile << " ";
      }
    }
    outfile << std::endl;
  }

  outfile.close();
  *verbose_out << "Saved confusion matrix to file " << filename << "." << std::endl;
}

void ForestClassification::writePredictionFile() {

  // Open prediction file for writing
  std::string filename = output_prefix + ".prediction";
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to prediction file: " + filename + ".");
  }

  // Write
  outfile << "Predictions: " << std::endl;
  if (predict_all) {
    for (size_t k = 0; k < num_trees; ++k) {
      outfile << "Tree " << k << ":" << std::endl;
      for (size_t i = 0; i < predictions.size(); ++i) {
        for (size_t j = 0; j < predictions[i].size(); ++j) {
          outfile << predictions[i][j][k] << std::endl;
        }
      }
      outfile << std::endl;
    }
  } else {
    for (size_t i = 0; i < predictions.size(); ++i) {
      for (size_t j = 0; j < predictions[i].size(); ++j) {
        for (size_t k = 0; k < predictions[i][j].size(); ++k) {
          outfile << predictions[i][j][k] << std::endl;
        }
      }
    }
  }

  *verbose_out << "Saved predictions to file " << filename << "." << std::endl;
}

void ForestClassification::saveToFileInternal(std::ofstream& outfile) {

  // Write num_variables
  outfile.write((char*) &num_variables, sizeof(num_variables));

  // Write treetype
  TreeType treetype = TREE_CLASSIFICATION;
  outfile.write((char*) &treetype, sizeof(treetype));

  // Write class_values
  saveVector1D(class_values, outfile);
}

void ForestClassification::loadFromFileInternal(std::ifstream& infile) {

  // Read number of variables
  size_t num_variables_saved;
  infile.read((char*) &num_variables_saved, sizeof(num_variables_saved));

  // Read treetype
  TreeType treetype;
  infile.read((char*) &treetype, sizeof(treetype));
  if (treetype != TREE_CLASSIFICATION) {
    throw std::runtime_error("Wrong treetype. Loaded file is not a classification forest.");
  }

  // Read class_values
  readVector1D(class_values, infile);

  for (size_t i = 0; i < num_trees; ++i) {

    // Read data
    std::vector<std::vector<size_t>> child_nodeIDs;
    readVector2D(child_nodeIDs, infile);
    std::vector<size_t> split_varIDs;
    readVector1D(split_varIDs, infile);
    std::vector<double> split_values;
    readVector1D(split_values, infile);

    // If dependent variable not in test data, change variable IDs accordingly
    if (num_variables_saved > num_variables) {
      for (auto& varID : split_varIDs) {
        if (varID >= dependent_varID) {
          --varID;
        }
      }
    }

    // Create tree
    Tree* tree = new TreeClassification(child_nodeIDs, split_varIDs, split_values, &class_values, &response_classIDs,
        &is_ordered_variable);
    trees.push_back(tree);
  }
}
// #nocov end
