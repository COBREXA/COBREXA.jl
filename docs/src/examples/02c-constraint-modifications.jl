
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

# # Making adjustments to the constraint system
#
# In the [previous example about model
# adjustments](02b-model-modifications.md), we noted that some constraint
# systems may be too complex to be changed within the limits of the usual FBC
# model view, and we may require a sharper tool to do the changes we need. This
# example shows how to do that by modifying the constraint systems that are
# generated within COBREXA to represent the metabolic model contents.

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import JSONFBCModels
import GLPK

model = load_model("e_coli_core.json") # flux balance type model

# ## Background: Constraint trees

# COBREXA uses [ConstraintTrees](https://github.com/COBREXA/ConstraintTrees.jl)
# to represent model structures internally. This framework provides a powerful
# unified interface over all constraints and variables in the model, making its
# manipulation much more convenient.

import ConstraintTrees as C

# In general, constraint-based models use fluxes as variables, and all the
# constraints are in terms of them (or derived quantities). We can get a
# constraint tree for the usual flux-balance-style models quite easily:

ct = flux_balance_constraints(model)

# The fluxes are represented by constraints for individual variables:

ct.fluxes

# The "mass balance" is represented as equality constraints:

ct.flux_stoichiometry

# The objective is represented as a "transparent reference" to the variables
# that specify the biomass. Notice that it has no bound (thus it's technically
# not a constraint, just a "label" for something that has a sensible semantic
# and can be constrained or optimized).

ct.objective

# ## Customizing the model

# You can easily create new values and constraints from the existing ones. For
# example, this is a total flux through exchanges of typical fermentation
# products:

total_fermentation = ct.fluxes.EX_ac_e.value + ct.fluxes.EX_etoh_e.value

# With the value in hand, we can constraint it (enforcing that the model
# outputs at least some fermentation products):

fermentation_constraint = C.Constraint(total_fermentation, (10.0, 1000.0))

# We can assign a name to the constraint, creating a small (singleton)
# constraint tree:

:fermentation^fermentation_constraint

# Named constraints can be freely combined, and we combine our new constraint
# with the whole original constraint tree, getting a differently constrained
# system:

fermenting_ct = ct * :fermentation^fermentation_constraint

# Constraint trees can be "solved", simply by choosing the objective and sending
# them to the appropriate function. Here, [`optimized_values`](@ref) rewrites
# the constraints into a JuMP model, which is subsequently solved and the
# solved variables are transformed back into semantically labeled values, in
# the same structure as the original constraint tree.

solution = optimized_values(
    fermenting_ct,
    objective = fermenting_ct.objective.value,
    optimizer = GLPK.Optimizer,
)

@test isapprox(solution.objective, 0.633738, atol = TEST_TOLERANCE) #src

# Models that can not be solved (for any reason) would instead return
# `nothing`. We demonstrate that by breaking the bounds of the original
# constraint trees to an unsolvable state:

ct.fluxes.ATPM.bound = C.Between(1000.0, 10000.0)

solution = optimized_values(ct, objective = ct.objective.value, optimizer = GLPK.Optimizer)

print(solution)

@test isnothing(solution) #src

# Several functions exist to simplify the construction of more complicated
# constraints. See the reference documentation for [generic constraint
# builders](reference/builders.md#Generic-constraints) for details.
