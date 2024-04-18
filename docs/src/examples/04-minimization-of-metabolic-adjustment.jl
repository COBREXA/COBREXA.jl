
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

# Minimization of metabolic adjustment analysis (MOMA) finds a flux solution
# that is closest, in an L2 sense, to some reference solution or model. Running
# this kind of analysis is straightforward using ConstraintTrees.

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import AbstractFBCModels.CanonicalModel as CM
import JSONFBCModels
import Clarabel

model = load_model("e_coli_core.json", CM.Model)

## Using reference solutions

# This is the simplest situation, where you supply a reference solution, and a
# flux solution that is closest to it is returned.

reference_fluxes =
    parsimonious_flux_balance_analysis(
        model,
        optimizer = Clarabel.Optimizer,
        settings = [silence],
    ).fluxes # get all the fluxes

modified_model = deepcopy(model) # change model, and avoid reference sharing

modified_model.reactions["CYTBD"]

modified_model.reactions["CYTBD"].upper_bound = 10.0

solution_with_reference_fluxes_l2 = metabolic_adjustment_minimization_analysis( # L2 norm used
    modified_model,
    reference_fluxes;
    optimizer = Clarabel.Optimizer,
    settings = [silence],
)

# It is also straightforward to use an L1 norm version of MOMA.

solution_with_reference_fluxes_l1 = linear_metabolic_adjustment_minimization_analysis(
    modified_model,
    reference_fluxes;
    optimizer = GLPK.Optimizer,
    settings = [silence],
)

@test isapprox(
    solution_with_reference_fluxes_l1.objective,
    0.21924874959,
    atol = TEST_TOLERANCE,
) #src
@test isapprox(solution_with_reference_fluxes_l1.fluxes.CYTBD, 10.0, atol = TEST_TOLERANCE) #src
@test solution_with_reference_fluxes_l1.linear_minimal_adjustment_objective < 305 #src


# ## Comparing models
# It is also possible to input metabolic models directly, saving some coding
# steps.

solution_l2 = metabolic_adjustment_minimization_analysis(
    modified_model,
    model;
    optimizer = Clarabel.Optimizer,
    settings = [silence],
)

@test isapprox(solution.objective, 0.241497, atol = TEST_TOLERANCE) #src
@test isapprox(
    solution_with_reference_fluxes_l2.objective,
    solution.objective,
    atol = TEST_TOLERANCE,
) #src
@test sqrt(solution.minimal_adjustment_objective) < 71 #src
@test isapprox(solution.fluxes.CYTBD, 10.0, atol = TEST_TOLERANCE) #src
@test isapprox( #src
    C.reduce( #src
        max, #src
        C.zip((a, b) -> abs(a - b), solution, solution_with_reference_fluxes, Float64), #src
        init = 0.0, #src
    ), #src
    0.0, #src
    atol = TEST_TOLERANCE, #src
) #src

# Likewise, the same thing can be done with the L1 norm MOMA

solution_l1 = linear_metabolic_adjustment_minimization_analysis(
    modified_model,
    model;
    optimizer = GLPK.Optimizer,
    settings = [silence],
)

@test isapprox(
    solution_l1.linear_minimal_adjustment_objective,
    304.409209299351,
    atol = TEST_TOLERANCE,
) #src

# ## Comparing fluxes

# Since we are using ConstraintTrees, it becomes trivial to compare logical
# units (e.g. groups of fluxes) between solutions. Here we will compare how
# different the fluxes are between the L1 and L2 norms. Since ConstraintTrees
# makes use of functional programming concepts, this is a simple one liner.

import ConstraintTrees as C

sort(
    collect(
        C.zip(
            -,
            solution_with_reference_fluxes_l2.fluxes,
            solution_with_reference_fluxes_l1.fluxes,
            Float64,
        ),
    ),
    by = last,
)
