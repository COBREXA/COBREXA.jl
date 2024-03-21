
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

# # Minimization of metabolic adjustment analysis

# TODO MOMA citation

using COBREXA

# TODO why not download+load combo?
download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import AbstractFBCModels.CanonicalModel as CM
import JSONFBCModels
import Clarabel

model = load_model("e_coli_core.json", CM.Model)

reference_fluxes =
    parsimonious_flux_balance_analysis(
        model,
        optimizer = Clarabel.Optimizer,
        settings = [silence],
    ).fluxes

model.reactions["CYTBD"]

model.reactions["CYTBD"].upper_bound = 10.0

solution = metabolic_adjustment_minimization_analysis(
    model,
    reference_fluxes;
    optimizer = Clarabel.Optimizer,
    settings = [silence],
)

solution.objective
@test isapprox(solution.objective, 0.241497, atol = TEST_TOLERANCE) #src

sqrt(solution.minimal_adjustment_objective)
@test sqrt(solution.minimal_adjustment_objective) < 71 #src

solution.fluxes.CYTBD
@test isapprox(solution.fluxes.CYTBD, 10.0, atol = TEST_TOLERANCE) #src

# compare

optimal_solution = parsimonious_flux_balance_analysis(
    model,
    optimizer = Clarabel.Optimizer,
    settings = [silence],
)

import ConstraintTrees as C
sort(collect(C.zip(-, optimal_solution.fluxes, solution.fluxes, Float64)), by = last)
