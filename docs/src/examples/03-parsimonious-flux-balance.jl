
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

# Here, we use [`parsimonious_flux_balance_analysis`](@ref) (pFBA) to find the
# optimal flux distribution in the *E. coli* "core" model. In essence, pFBA
# first uses FBA to find an optimal objective value for the model, and then
# minimizes the squared distance of the flux from the zero (i.e., minimizes its
# L2 norm). As the main benefit, this gives a unique (and possibly more
# realistic) solution to the model.

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

# Notably, we need an optimizer that can solve quadratic (QP) models:
import Clarabel

import JSONFBCModels
model = load_model("e_coli_core.json")

# The pFBA is, in its most default form, implemented in function
# [`parsimonious_flux_balance_analysis`](@ref):

solution = parsimonious_flux_balance_analysis(
    model;
    optimizer = Clarabel.Optimizer,
    settings = [silence],
)

#

solution.fluxes

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

# ## Using different solvers for the problem stages
#
# Sometimes it is useful to employ a dedicated LP solver to find the solution to
# the original FBA problem, and then a dedicated QP solver to minimize the
# fluxes. We can set the optimizer and parsimonious optimizer separately using
# keyword arguments:

import GLPK

solution = parsimonious_flux_balance_analysis(
    model;
    optimizer = GLPK.Optimizer, # GLPK is good for LP but cannot do QP
    settings = [silence],
    parsimonious_optimizer = Clarabel.Optimizer, # Clarabel is not very precise but can solve QP
)

@test isapprox(solution.objective, 0.873922; atol = TEST_TOLERANCE) #src

# ## Using linear parsimonious

# For efficiency reasons, it is also possible to use a pFBA version that
# optimizes the L1 norm instead of the L2 one (i.e., minimizing a sum of
# absolute values instead of the sum of squares). In turn, the uniqueness
# property of the solution is lost, but we do not need a QP-capable optimizer
# at all:

linear_solution =
    linear_parsimonious_flux_balance_analysis(model; optimizer = GLPK.Optimizer)

#

linear_solution.fluxes

@test isapprox(linear_solution.parsimonious_objective, 518.422; atol = TEST_TOLERANCE) #src
@test isapprox( #src
    sum(abs.(values(linear_solution.fluxes))), #src
    linear_solution.parsimonious_objective, #src
    atol = TEST_TOLERANCE, #src
) #src
