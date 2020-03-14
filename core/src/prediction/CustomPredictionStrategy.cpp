/*-------------------------------------------------------------------------------
  This file is part of generalized random forest (grf).

  grf is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  grf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with grf. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include <cmath>
#include <string>
#include "prediction/CustomPredictionStrategy.h"

const size_t CustomPredictionStrategy::OUTCOME = 0;             // use to record the numerator of estimator
const std::size_t CustomPredictionStrategy::TREATMENT = 1;      // use to record the denumerator of estimator
const std::size_t CustomPredictionStrategy::INSTRUMENT = 2;

const std::size_t NUM_TYPES = 2;


size_t CustomPredictionStrategy::prediction_length() {
    return 1;
}

// The estimator has nominator part and denominator part !
// See OptimizedPredictionCollector line143
std::vector<double> CustomPredictionStrategy::predict(const std::vector<double>& average) {
    return { average.at(OUTCOME) / average.at(TREATMENT) };
}


// leaf_values contains leaf_values(the nominator, denominator of each trees), num_trees, strategy->prediction_value_length()
/*
std::vector<double> CustomPredictionStrategy::compute_variance(
                                                                   const std::vector<double>& average,
                                                                   const PredictionValues& leaf_values,
                                                                   uint ci_group_size) {
    
    double average_outcome = average.at(OUTCOME) / average.at(TREATMENT);
    
    double num_good_groups = 0;
    double psi_squared = 0;
    double psi_grouped_squared = 0;
    
    // leaf_values.get_num_nodes() == num_trees
    for (size_t group = 0; group < leaf_values.get_num_nodes() / ci_group_size; ++group) {
        bool good_group = true;
        for (size_t j = 0; j < ci_group_size; ++j) {
            if (leaf_values.empty(group * ci_group_size + j)) {
                good_group = false;
            }
        }
        if (!good_group) continue;
        
        num_good_groups++;
        
        double group_psi = 0;
        
        for (size_t j = 0; j < ci_group_size; ++j) {
            size_t i = group * ci_group_size + j;
            double psi_1 = leaf_values.get(i, OUTCOME) / leaf_values.get(i, TREATMENT) - average_outcome;
            
            psi_squared += psi_1 * psi_1;
            group_psi += psi_1;
        }
        
        group_psi /= ci_group_size;
        psi_grouped_squared += group_psi * group_psi;
    }
    
    double var_between = psi_grouped_squared / num_good_groups;
    double var_total = psi_squared / (num_good_groups * ci_group_size);
    
    // This is the amount by which var_between is inflated due to using small groups
    double group_noise = (var_total - var_between) / (ci_group_size - 1);
    
    // A simple variance correction, would be to use:
    // var_debiased = var_between - group_noise.
    // However, this may be biased in small samples; we do an objective
    // Bayes analysis of variance instead to avoid negative values.
    double var_debiased = bayes_debiaser.debias(var_between, group_noise, num_good_groups);
    
    
    //leaf_values.get_values(0)
    //std::vector<double> vect(7, -1);
    return { var_debiased };

}
*/

size_t CustomPredictionStrategy::prediction_value_length() {
    return NUM_TYPES;
}



// prediction value for each leaf of a tree
// For each leaf, averages[OUTCOME] = sum of IPCW *Z; averages[TREATMENT] = sum of IPCW; averages[INSTRUMENT] = number of data points
PredictionValues CustomPredictionStrategy::precompute_prediction_values(
                                                                        const std::vector<std::vector<size_t>>& leaf_samples,
                                                                        const Observations& observations) {
    size_t num_leaves = leaf_samples.size();
    std::vector<std::vector<double>> values(num_leaves); //values[i] is values of nodes in leaf i
    
    for (size_t i = 0; i < num_leaves; i++) {
        const std::vector<size_t>& leaf_node = leaf_samples.at(i);
        if (leaf_node.empty()) {
            continue;
        }
        
        std::vector<double>& averages = values[i];  // Note that if leaf_node is empty, then values[i] will be empty !!!
        averages.resize(NUM_TYPES);
        
        double muhat = 0.0;
        double sum_weights=0.0;
        double sample_counter=0.0;
        
        for (auto& sample : leaf_node) {
            double Z= observations.get(Observations::OUTCOME, sample);
            double delta= observations.get(Observations::TREATMENT, sample);
            double G= observations.get(Observations::INSTRUMENT, sample);
            double IPCW = delta/(1-G);
            muhat += IPCW*Z;
            sum_weights+=IPCW;
            sample_counter+=1;
        }
        //averages[INSTRUMENT] = sample_counter;
        if(sum_weights>0){
            averages[OUTCOME] = muhat/sample_counter;
            averages[TREATMENT] = sum_weights/sample_counter;
        }else{
            averages[OUTCOME] = 0;
            averages[TREATMENT] = 0;
        }
        
    }
    return PredictionValues(values, num_leaves, NUM_TYPES);
}

std::vector<double> CustomPredictionStrategy::compute_debiased_error(
                                                                         size_t sample,
                                                                         const std::vector<double>& average,
                                                                         const PredictionValues& leaf_values,
                                                                         const Observations& observations) {
    /*
    double outcome = observations.get(Observations::OUTCOME, sample);
    
    double error = average.at(OUTCOME) - outcome;
    double mse = error * error;
    
    double bias = 0.0;
    size_t num_trees = 0;
    for (size_t n = 0; n < leaf_values.get_num_nodes(); n++) {
        if (leaf_values.empty(n)) {
            continue;
        }
        
        double tree_variance = leaf_values.get(n, OUTCOME) - average.at(OUTCOME);
        bias += tree_variance * tree_variance;
        num_trees++;
    }
    
    if (num_trees <= 1) {
        return { NAN };
    }
    
    bias /= num_trees * (num_trees - 1);
    return { mse - bias };
     */
    return {0};
}









/*
    PredictionValues CustomPredictionStrategy::precompute_prediction_values(
                                                                            const std::vector<std::vector<size_t>>& leaf_samples,
                                                                            const Observations& observations) {
        size_t num_leaves = leaf_samples.size();
        std::vector<std::vector<double>> values(num_leaves); //values[i] is values of nodes in leaf i
        
        for (size_t i = 0; i < num_leaves; i++) {
            const std::vector<size_t>& leaf_node = leaf_samples.at(i);
            if (leaf_node.empty()) {
                continue;
            }
            
            std::vector<double>& averages = values[i];
            averages.resize(1);
            
            double muhat = 0.0;
            double sum_weights=0.0;
            double check=0.0;
            int count=-3;
            
            for (auto& sample : leaf_node) {
                double Z= observations.get(Observations::OUTCOME, sample);
                double delta= observations.get(Observations::TREATMENT, sample);
                double G= observations.get(Observations::INSTRUMENT, sample);
                double IPCW = delta/(1-G);
                muhat += IPCW*Z;
                sum_weights+=IPCW;
            }
            
            if(sum_weights>0){
                averages[OUTCOME] = muhat/sum_weights;
            }else{
                averages[OUTCOME] = 0;
            }
            
        }
        return PredictionValues(values, num_leaves, 1);
    }
*/


std::vector<double> CustomPredictionStrategy::compute_variance(
                                                               const std::vector<double>& average,
                                                               const PredictionValues& leaf_values,
                                                               uint ci_group_size) {
    
    
    double average_outcome = average.at(OUTCOME)/average.at(TREATMENT);
    
    double num_good_groups = 0;
    double psi_squared = 0;
    double psi_grouped_squared = 0;
    
    // leaf_values.get_num_nodes() == num_trees
    for (size_t group = 0; group < leaf_values.get_num_nodes() / ci_group_size; ++group) {
        bool good_group = true;
        for (size_t j = 0; j < ci_group_size; ++j) {
            if (leaf_values.empty(group * ci_group_size + j)) {
                good_group = false;
            }
        }
        if (!good_group) continue;
        
        num_good_groups++;
        
        double group_psi = 0;
        
        for (size_t j = 0; j < ci_group_size; ++j) {
            size_t i = group * ci_group_size + j;
            double psi_1 = leaf_values.get(i, OUTCOME)-average_outcome*leaf_values.get(i, TREATMENT);
            psi_squared += psi_1 * psi_1;
            group_psi += psi_1;
        }
        
        group_psi /= ci_group_size;
        psi_grouped_squared += group_psi * group_psi;
    }
    
    double var_between = psi_grouped_squared / num_good_groups;
    double var_total = psi_squared / (num_good_groups * ci_group_size);
    
    // This is the amount by which var_between is inflated due to using small groups
    double group_noise = (var_total - var_between) / (ci_group_size - 1);
    
    // A simple variance correction, would be to use:
    // var_debiased = var_between - group_noise.
    // However, this may be biased in small samples; we do an objective
    // Bayes analysis of variance instead to avoid negative values.
    double var_debiased = bayes_debiaser.debias(var_between, group_noise, num_good_groups);
    
    
    //leaf_values.get_values(0)
    //std::vector<double> vect(7, -1);
    return { var_debiased };
}
