
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

# # Parsimonious flux balance analysis

# We will use [`parsimonious_flux_balance_analysis`](@ref) to find the optimal
# flux distribution in the *E. coli* "core" model.
#
# TODO pFBA citation

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

# next, load the necessary packages

import JSONFBCModels
import Clarabel # can solve QPs

model = load_model("e_coli_core.json") # load the model

# Use the convenience function to run standard pFBA on

solution = parsimonious_flux_balance_analysis(
    model;
    optimizer = Clarabel.Optimizer,
    settings = [silence],
)

@test isapprox(solution.objective, 0.873922; atol = TEST_TOLERANCE) #src
@test isapprox( #src
    solution.parsimonious_objective, #src
    11414.211988071253, #src
    atol = QP_TEST_TOLERANCE, #src
) #src
@test isapprox( #src
    sum(x^2 for x in values(solution.fluxes)), #src
    solution.parsimonious_objective, #src
    atol = QP_TEST_TOLERANCE, #src
) #src

# ## Linear version
#
# TODO

import GLPK

linear_solution =
    linear_parsimonious_flux_balance_analysis(model; optimizer = GLPK.Optimizer)

@test isapprox(linear_solution.objective, 0.873922; atol = TEST_TOLERANCE) #src
@test isapprox(linear_solution.parsimonious_objective, 518.422; atol = TEST_TOLERANCE) #src
@test isapprox( #src
    sum(abs.(values(linear_solution.fluxes))), #src
    linear_solution.parsimonious_objective, #src
    atol = TEST_TOLERANCE, #src
) #src
