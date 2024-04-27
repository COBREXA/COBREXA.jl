
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

# # Flux variability analysis (FVA)

# FVA finds a range of fluxes through each reaction where the model can behave
# optimally. In brief, it first runs a FBA to get the optimal objective value,
# constraints the model to the optimal (or near-optimal) space, and runs a
# separate minimization and maximization task for each reaction to find their
# individual ranges.

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import JSONFBCModels, GLPK

model = load_model("e_coli_core.json")

# The "usual" form of FBA is available via the eponymous function:

solution = flux_variability_analysis(model, optimizer = GLPK.Optimizer)

@test isapprox(solution.ACALD[1], -2.542370370370188, atol = TEST_TOLERANCE) #src
@test isapprox(solution.ACALD[2], 0.0, atol = TEST_TOLERANCE) #src

# ## Specifying objective bounds

# By default, FVA computes variability from the feasible region that is bounded
# to be within 10% of the optimal objective value. That is not very strict, and
# you can force much lower tolerance.

# Here, we force the optimal region to be within 0.00001 units of the optimal
# objective value:

very_close = flux_variability_analysis(
    model,
    optimizer = GLPK.Optimizer,
    objective_bound = absolute_tolerance_bound(1e-5),
)

# Here, we relax that to 1% of the optimal objective value:

one_percent_close = flux_variability_analysis(
    model,
    optimizer = GLPK.Optimizer,
    objective_bound = relative_tolerance_bound(0.99),
)

#md # !!! tip "Speed up FVA with parallel processing"
#md #     By default, FVA is parallelized on all workers that are available in the worker pool of the `Distributed` package, which may speed up the computation considerably. See the [parallel processing documentation](../distributed.md) for more details.
