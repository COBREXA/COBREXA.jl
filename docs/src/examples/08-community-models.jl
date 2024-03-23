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

# # Community FBA models

using COBREXA

# Here we will construct a community FBA model of two  *E. coli* "core" models
# that can interact by exchanging selected metabolites. To do this, we will need
# the model, which we can download if it is not already present.

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

# Additionally to COBREXA and the model format package, we will need a solver
# -- let's use GLPK here:

import JSONFBCModels
import GLPK
import AbstractFBCModels.CanonicalModel as CM

ecoli1 = load_model("e_coli_core.json", CM.Model)
ecoli1.reactions["EX_glc__D_e"].lower_bound = -1000.0
ecoli1.reactions["EX_glc__D_e"].upper_bound = 1000.0
ecoli2 = deepcopy(ecoli1)

# customize models a bit

ecoli1.reactions["CYTBD"].lower_bound = ecoli1.reactions["CYTBD"].upper_bound = 0.0
ecoli2.reactions["FBA"].lower_bound = ecoli2.reactions["FBA"].upper_bound = 0.0

# solve

solution = community_flux_balance_analysis(
    ["bug1" => (ecoli1, 0.2), "bug2" => (ecoli2, 0.8)],
    ["EX_glc__D_e" => (-10.0, 0.0)],
    optimizer = GLPK.Optimizer,
)

@test isapprox(solution.community_biomass, 0.5237157737585179, atol = TEST_TOLERANCE) #src

@test isapprox( #src
    solution.bug1.fluxes.BIOMASS_Ecoli_core_w_GAM, #src
    solution.bug2.fluxes.BIOMASS_Ecoli_core_w_GAM, #src
    atol = TEST_TOLERANCE, #src
) #src

# ## Compare

# We can now e.g. observe the differences in individual pairs of exchanges:

import ConstraintTrees as C

C.zip(
    tuple,
    solution.bug1.interface.exchanges,
    solution.bug2.interface.exchanges,
    Tuple{Float64,Float64},
)

# ## Which composition is best?

screen(0.0:0.1:1.0) do ratio2
    ratio1 = 1 - ratio2
    res = community_flux_balance_analysis(
        ["bug1" => (ecoli1, ratio1), "bug2" => (ecoli2, ratio2)],
        ["EX_glc__D_e" => (-10.0, 0.0)],
        optimizer = GLPK.Optimizer,
    )
    (ratio1, ratio2) => (isnothing(res) ? nothing : res.community_biomass)
end

# ...seems a lot like `bug1` will eventually disappear.

# # TODO old version

# Community models work by joining its members together through their exchange
# reactions, weighted by the abundance of each microbe. These exchange reactions
# are then linked to an environmental exchange. For more theoretical details,
# see "Gottstein, et al, 2016, Constraint-based stoichiometric modelling from
# single organisms to microbial communities, Journal of the Royal Society
# Interface".

# ## Building a community of two *E. coli*s

# Here we will construct a simple community of two interacting microbes. To do
# this, we need to import the models. We will represent the models only as
# constraint trees, because it is easier to build the model explicitly than
# rely on an opaque one-shot function.

# ecoli1 = flux_balance_constraints(model, interface = :sbo)
#ecoli2 = flux_balance_constraints(model, interface = :sbo)

# Since the models are usually used in a mono-culture context, the glucose input
# for each individual member is limited. We need to undo this limitation, and
# rather rely on the constrained environmental exchange reaction (and the bounds
# we set for it earlier).
#ecoli1.fluxes.EX_glc__D_e.bound = C.Between(-1000.0, 1000.0)
#ecoli2.fluxes.EX_glc__D_e.bound = C.Between(-1000.0, 1000.0)

# To make the community interesting, we can limit different reactions in both
# members to see how the  models cope together:
#ecoli1.fluxes.CYTBD.bound = C.Between(-10.0, 10.0)
#ecoli2.fluxes.ACALD.bound = C.Between(-5.0, 5.0)

# Because we created the trees with interfaces, we can connect them easily to
# form a new model with the interface. For simplicity, we use the
# interface-scaling functionality of [`interface_constraints`](@ref
# COBREXA.interface_constraints) to bring in cFBA-like community member
# abundances:

#cc = interface_constraints(
#:bug1 => (ecoli1, ecoli1.interface, 0.2),
#:bug2 => (ecoli2, ecoli2.interface, 0.8),
#)

# To make the community behave as expected, we need to force equal (scaled)
# growth of all members:

#cc *=
#:equal_growth^equal_value_constraint(
#cc.bug1.fluxes.BIOMASS_Ecoli_core_w_GAM,
#cc.bug2.fluxes.BIOMASS_Ecoli_core_w_GAM,
#)
#
# Now we can simulate the community growth by optimizing the new "interfaced"
# biomass:

#optimized_cc = optimized_values(
#cc,
#objective = cc.interface.biomass.BIOMASS_Ecoli_core_w_GAM.value,
#optimizer = GLPK.Optimizer,
#)
