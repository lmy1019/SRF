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

#include "CustomRelabelingStrategy.h"

std::unordered_map<size_t, double> CustomRelabelingStrategy::relabel(
    const std::vector<size_t>& samples,
    const Observations& observations) {
    
    std::unordered_map<size_t, double> relabeled_observations;
    double sum_IPCW=0;
    double muhat=0;
    double node_size = 0;
    
    for (size_t sample : samples) {
        double outcome = observations.get(Observations::OUTCOME, sample);
        double delta = observations.get(Observations::TREATMENT, sample);
        double G = observations.get(Observations::INSTRUMENT, sample);
        double IPCW = 1.0*delta/(1-G);
        muhat += outcome*IPCW;
        sum_IPCW += IPCW;
        node_size+=1;
    }
    muhat = muhat/sum_IPCW;
    
    for (size_t sample : samples) {
        double outcome = observations.get(Observations::OUTCOME, sample);
        double delta = observations.get(Observations::TREATMENT, sample);
        double G = observations.get(Observations::INSTRUMENT, sample);
        double IPCW = 1.0*delta/(1-G);
        relabeled_observations[sample]  = IPCW*(outcome-muhat)/(sum_IPCW/node_size);
    }
    return relabeled_observations;
}
