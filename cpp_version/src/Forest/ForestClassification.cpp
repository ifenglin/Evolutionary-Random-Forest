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
	  if (!class_count.empty()) {
		  predictions[0][0][sample_idx] = mostFrequentValue(class_count, random_number_generator);
	  }
	  else {
		  predictions[0][0][sample_idx] = 0;
	  }
    }
  }
}

void ForestClassification::computePredictionErrorInternal() {

	// Use a fraction of the samples
	//size_t num_samples_downsized = (size_t)sqrt(num_samples);
	size_t num_samples_downsized = (size_t)(num_samples);
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
  *verbose_out << "--Recording predictions.." << std::endl;
  //predictions_each_tree = std::vector<std::vector<std::vector<double>>>(1, std::vector<std::vector<double>>(num_trees, std::vector<double>(num_samples)));
  std::vector<size_t> num_missclassifications_each_tree = std::vector<size_t>(num_trees, 0);
  std::vector<size_t> num_predictions_each_tree = std::vector<size_t>(num_trees, 0);
  //prediction_error_each_tree = std::vector<double>(num_trees);
  try {
	  for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
		  for (size_t sample_idx = 0; sample_idx < num_samples_downsized; ++sample_idx) {
			  size_t sampleID = sampleIDs[sample_idx];
			  size_t sample_idx_tree = trees[tree_idx]->findOobSample(sampleID);
			  // if the sample is in the tree's Obb samples
			  if (sample_idx_tree != trees[tree_idx]->getNumSamplesOob()) {
				  double value = ((TreeClassification*)trees[tree_idx])->getPrediction(sample_idx_tree);
				  double real_value = data->get(sample_idx, dependent_varID);
				  /*if (class_counts[sampleID].find(value) == class_counts[sampleID].end()) {
					  class_counts[sampleID][value] = 0;
				  }*/
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
  *verbose_out << "--Computing majority.." << std::endl;
  predictions = std::vector<std::vector<std::unordered_map<size_t, unsigned>>>(1, std::vector<std::unordered_map<size_t, unsigned>>(1, std::unordered_map<size_t, unsigned>(num_samples_downsized)));
  std::unordered_map<size_t, double> secondPredictions;
  try {
	  for (size_t i = 0; i < num_samples_downsized; ++i) {
		  size_t sampleID = sampleIDs[i];
		  if (!class_counts[sampleID].empty()) {
			  predictions[0][0][sampleID] = mostFrequentValue(class_counts[sampleID], random_number_generator);
			  if (class_counts[sampleID].size() > 1) {
				  double real_value = data->get(sampleID, dependent_varID);
				  secondPredictions[sampleID] = mostFrequentFalseValue(class_counts[sampleID], real_value);
			  }
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

  // Compare forest predictions with true values
  *verbose_out << "--Comparing predictions.." << std::endl;
  size_t num_missclassifications = 0;
  size_t num_predictions = 0;
  double sum_margin = 0;
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
			  sum_margin += computeMargin(class_counts[sampleID], real_value);
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
  overall_strength = std::max( 0.0, sum_margin / (double) num_predictions);
  overall_correlation = 0;
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
  *verbose_out << "--Calculating correlation.." << std::endl;
  try {
	  correlation_each_tree = std::vector<std::vector<std::vector<double>>>(1, std::vector<std::vector<double>>(num_trees, std::vector<double>(num_trees, 0.0)));
	  //std::vector<size_t> intersected_oob_sampleIDs;
	  //std::vector<size_t> i_oob_sampleIDs;
	  //std::vector<size_t> j_oob_sampleIDs;
	  for ( size_t i = 0; i < num_trees; ++i) {
		  //i_oob_sampleIDs = trees[i]->getOobSampleIDs();
		  for ( size_t j = i + 1; j < num_trees; ++j) {
			size_t correlated_count = 0;
			size_t uncorrelated_count = 0;
			double n = 0, sum_xy = 0, sum_x = 0, sum_y = 0, sum_x2 = 0, sum_y2 = 0, sum_true_true = 0;
			//j_oob_sampleIDs = trees[j]->getOobSampleIDs();
			//intersected_oob_sampleIDs.clear();
			//std::set_intersection(i_oob_sampleIDs.begin(), i_oob_sampleIDs.end(), j_oob_sampleIDs.begin(), j_oob_sampleIDs.end(), std::back_inserter(intersected_oob_sampleIDs));
			for ( size_t sample_idx = 0; sample_idx < num_samples_downsized; ++sample_idx) {
				size_t sampleID = sampleIDs[sample_idx];
				// check if sampleID is in oob sample in both trees
				// stage-wise checking reduces the times of searching
				// equivalent to the following, but assumably faster when oob samples are more than non-oob samples.
				// if(std::find(intersected_oob_sampleIDs.begin(), intersected_oob_sampleIDs.end(), sampleID) != intersected_oob_sampleIDs.end()){}
				size_t sample_idx_tree_i = trees[i]->findOobSample(sampleID);
				if (sample_idx_tree_i != trees[i]->getNumSamplesOob()) {
					size_t sample_idx_tree_j = trees[j]->findOobSample(sampleID);
					if (sample_idx_tree_j != trees[j]->getNumSamplesOob()) {
						//if (std::find(intersected_oob_sampleIDs.begin(), intersected_oob_sampleIDs.end(), sampleID) != intersected_oob_sampleIDs.end()) {
						double true_value = data->get(sampleID, dependent_varID);
						double predicted_value_i = ((TreeClassification*)trees[i])->getPrediction(sample_idx_tree_i);
						double predicted_value_j = ((TreeClassification*)trees[j])->getPrediction(sample_idx_tree_j);
						double x, y, best_false_values;

						if (secondPredictions.find(sampleID) != secondPredictions.end()) {
							best_false_values = secondPredictions[sampleID];
							x = (int)(predicted_value_i == true_value) - (int)(predicted_value_i == best_false_values);
							y = (int)(predicted_value_j == true_value) - (int)(predicted_value_j == best_false_values);

							// skip if both predictions are true class.
							if (x == 1 || y == 1) {
								++sum_true_true;
							}
							else {
								sum_xy += x * y;
								sum_x += x;
								sum_y += y;
								sum_x2 += abs(x);
								sum_y2 += abs(y);
								n++;
							}
						}
						
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
			}

			size_t total_count = correlated_count + uncorrelated_count;
			double correlation_rate = 0.0, no_true_correlation_rate = 0.0;
			if (n != 0) {
				// calculate correlations without true-true samples
				// correlation_rate = (double)correlated_count / (double)(total_count);
				
				no_true_correlation_rate = ((n * sum_xy) - (sum_x * sum_y)) /
					(sqrt((n * sum_x2) - pow(sum_x, 2)) * sqrt((n * sum_y2) - pow(sum_y, 2)));

				if ( std::isinf(no_true_correlation_rate) ) {
					// correlation is underfined when one of the variances is zero
					// as zero variance is unwanted for a tree, use 1 for its value;
					no_true_correlation_rate = 1.0;
				}

				// transform correlation into [0, 1] range.
				no_true_correlation_rate = (no_true_correlation_rate + 1.0) / 2.0;

				// add true predictions into sum for Breiman's correlation
				sum_xy += sum_true_true;
				sum_x += sum_true_true;
				sum_y += sum_true_true;
				sum_x2 += sum_true_true;
				sum_y2 += sum_true_true;
				n += sum_true_true;
				// calculate correlations with true-ture samples;
				correlation_rate = ((n * sum_xy) - (sum_x * sum_y)) /
					(sqrt((n * sum_x2) - pow(sum_x, 2)) * sqrt((n * sum_y2) - pow(sum_y, 2)));
				
				if ( std::isinf(correlation_rate) ) {
					// correlation is underfined when one of the variances is zero
					// as zero variance is unwanted for a tree, use 1 for its value;
					correlation_rate = 1.0;
				}
			}

			// set correlation only to the tree that has a higher prediction error
			if (prediction_error_each_tree[i] < prediction_error_each_tree[j]) {
				correlation_each_tree[0][i][j] = 0.0;
				correlation_each_tree[0][j][i] = no_true_correlation_rate;
			}
			else if (prediction_error_each_tree[i] > prediction_error_each_tree[j]) {
				correlation_each_tree[0][i][j] = no_true_correlation_rate;
				correlation_each_tree[0][j][i] = 0.0;
			}
			else {
				correlation_each_tree[0][i][j] = no_true_correlation_rate;
				correlation_each_tree[0][j][i] = no_true_correlation_rate;
			}
			overall_correlation += correlation_rate;
		  }
	  }
	  overall_correlation /= num_trees * (num_trees + 1) / 2;
  } catch (const std::bad_alloc &e) {
	  std::cout << "Allocation failed when calculating correation between tree: " << e.what() << ". Skipped" << std::endl;
  }

  margin_each_tree = std::vector<double>(num_trees, 0);
  // calculate variance
  *verbose_out << "--Calculating variance.." << std::endl;
  try {
	  for (size_t i = 0; i < num_trees; ++i) {
	      double n = 0, sum_x = 0, sum_x2 = 0;
		  for (size_t sample_idx = 0; sample_idx < num_samples_downsized; ++sample_idx) {
			  size_t sampleID = sampleIDs[sample_idx];
			  size_t sample_idx_tree = trees[i]->findOobSample(sampleID);
	     	  if (sample_idx_tree != trees[i]->getNumSamplesOob()) {
				  double true_value = data->get(sampleID, dependent_varID);
				  double predicted_value = ((TreeClassification*)trees[i])->getPrediction(sample_idx_tree);
				  double x,  best_false_values;
				  if (secondPredictions.find(sampleID) != secondPredictions.end()) {
				      best_false_values = secondPredictions[sampleID];
				      x = (int)(predicted_value == true_value) - (int)(predicted_value == best_false_values);
				  	  sum_x += x;
					  sum_x2 += abs(x);
					  n++;
				  }
			  } 
		  }
		  size_t variance = (n * sum_x2) - pow(sum_x, 2);
		  overall_variance += variance;
		  //margin_each_tree[i] = sum_x / n;
		  // define new margin as the inverse of variance in 1e6
		  margin_each_tree[i] = ( (sum_x > 0)? 1.0 : -1.0 ) / variance * 1000000.0;
	  }
	  overall_variance /= num_trees;
  }
  catch (const std::bad_alloc &e) {
	  std::cout << "Allocation failed when calculating correation between tree: " << e.what() << ". Skipped" << std::endl;
  }
}

// #nocov start
void ForestClassification::writeOutputInternal() {
  *verbose_out << "Tree type:                         " << "Classification" << std::endl;
}

void ForestClassification::writeConfusionFile() {

  // Open confusion file for writing
  std::string filename = output_prefix + "_confusion.txt";
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
  std::string filename = output_prefix + "_prediction.txt";
  std::ofstream outfile;
  if (prediction_overwrite) {
	  outfile.open(filename, std::ios::out);
	  prediction_overwrite = false;
  }
  else {
	  outfile.open(filename, std::ios::app);
  }
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to prediction file: " + filename + ".");
  }

  // Write
  //outfile << "Predictions: " << std::endl;
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
