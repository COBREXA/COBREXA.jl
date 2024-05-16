
# Copyright (c) 2021-2024, University of Luxembourg                         #src
# Copyright (c) 2021-2024, Heinrich-Heine University Duesseldorf            #src
#                                                                           #src
# Licensed under the Apache License, Version 2.0 (the "License");           #src
# you may not use this file except in compliance with the License.          #src
# You may obtain a copy of the License at                                   #src
#                                                                           #src
#     http://www.apache.org/licenses/LICENSE-2.0                            #src
#                                                                           #src
# Unless required by applicable law or agreed to in writing, software       #src
# distributed under the License is distributed on an "AS IS" BASIS,         #src
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  #src
# See the License for the specific language governing permissions and       #src
# limitations under the License.                                            #src

# # Flux sampling
#
# Flux sampling gives an interesting statistical insight into the behavior of
# the model in the optimal feasible space, and the general "shape" of the
# optimal- or near-optimal set of feasible states of a given model.

# For demonstration, we need the usual packages and models:

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import JSONFBCModels, GLPK

model = load_model("e_coli_core.json")

# Function [`flux_sample`](@ref) uses linear optimization to generate a set of
# warm-up points (by default, the method to generate the warm-up is basically
# FVA), and then runs the hit-and-run flux sampling algorithm on the
# near-optimal feasible space of the model:
s = flux_sample(
    model,
    optimizer = GLPK.Optimizer,
    objective_bound = relative_tolerance_bound(0.99),
    n_chains = 2,
    collect_iterations = [10],
)

@test 21.8 < sum(s.O2t) / length(s.O2t) < 22.0 #src

# The result is a tree of vectors of sampled states for each value; the order
# of the values in these vectors is fixed. You can thus e.g. create a good
# matrix for plotting the sample as 2D scatterplot:

[ s.O2t s.CO2t ]
