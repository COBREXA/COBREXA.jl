
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

# FVA performs FBA to find an optimal objective, and then constrains this
# objective to lie within certain bounds. Thereafter, each other reaction is
# minimized and maximized to identify its maximum variability.

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import JSONFBCModels, GLPK

model = load_model("e_coli_core.json")

# Most of the basic analysis functions have standardized frontend functions,
# making it easy to call them.

solution = flux_variability_analysis(model, optimizer = GLPK.Optimizer)

@test isapprox(solution.ACALD[1], -2.542370370370188, atol = TEST_TOLERANCE) #src
@test isapprox(solution.ACALD[2], 0.0, atol = TEST_TOLERANCE) #src

# ## Specifying bounds

# Options for FVA include changing the variability tolerances on the objective,
# or making use of parallel processing, to speed up computations.

very_close = flux_variability_analysis(
    model,
    optimizer = GLPK.Optimizer,
    objective_bound = absolute_tolerance_bound(1e-5),
)

one_percent_close = flux_variability_analysis(
    model,
    optimizer = GLPK.Optimizer,
    objective_bound = relative_tolerance_bound(0.99),
)

using Distributed
addprocs(2) # add workers to distribute optimization problem across more CPUs
@everywhere using COBREXA, GLPK # load packages on workers

solution = flux_variability_analysis(model, optimizer = GLPK.Optimizer; workers=workers()) # distribute work to more cores
