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

# Here we construct a community FBA model of two  *E. coli* "core" models that
# can interact by exchanging selected metabolites. To do this, we will need the
# model, which we can download if it is not already present.

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

# Additionally to COBREXA and the model format package, we will need a solver
# and a few supporting packages.

import JSONFBCModels
import HiGHS
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C

# The core model has an artificial bound on input glucose; here we unblock that
# one, and we are going to add a community-global glucose intake bound later.

ecoli = load_model("e_coli_core.json", CM.Model)
ecoli.reactions["EX_glc__D_e"].lower_bound = -1000.0
ecoli.reactions["EX_glc__D_e"].upper_bound = 1000.0

# To create a community that is actually interesting, we need some diversity.
# Here we simply block a different reaction in each of the community members:

ecoli1 = deepcopy(ecoli)
ecoli1.reactions["CYTBD"].lower_bound = ecoli1.reactions["CYTBD"].upper_bound = 0.0
ecoli2 = deepcopy(ecoli)
ecoli2.reactions["FBA"].lower_bound = ecoli2.reactions["FBA"].upper_bound = 0.0

# ## Analysing the community

# To construct the community, we have to provide identifiers for the models
# (these will be used in the constraint tree), and corresponding models with
# the abundances.
my_community = Dict("bug1" => (ecoli1, 0.2), "bug2" => (ecoli2, 0.8))

# The community is constructed and analysed using
# [`community_flux_balance_analysis`](@ref):
solution = community_flux_balance_analysis(
    my_community,
    ["EX_glc__D_e" => (-10.0, 0.0)],
    optimizer = HiGHS.Optimizer,
)

@test isapprox(solution.community_biomass, 0.5237157737585179, atol = TEST_TOLERANCE) #src

@test isapprox( #src
    solution.bug1.fluxes.BIOMASS_Ecoli_core_w_GAM, #src
    solution.bug2.fluxes.BIOMASS_Ecoli_core_w_GAM, #src
    atol = TEST_TOLERANCE, #src
) #src

# ## Investigating the solution

# We can now e.g. observe the differences in individual pairs of exchanges:

C.zip(
    tuple,
    solution.bug1.interface.exchanges,
    solution.bug2.interface.exchanges,
    Tuple{Float64,Float64},
)

# ...or use [`screen`](@ref) to efficiently find out which composition is best:

screen(0.0:0.1:1.0) do ratio2
    ratio1 = 1 - ratio2
    res = community_flux_balance_analysis(
        [("bug1" => (ecoli1, ratio1)), ("bug2" => (ecoli2, ratio2))],
        ["EX_glc__D_e" => (-10.0, 0.0)],
        interface = :sbo, # usually more reproducible
        optimizer = HiGHS.Optimizer,
    )
    (ratio1, ratio2) => (isnothing(res) ? nothing : res.community_biomass)
end

# (The result seem like the `bug1` is eventually going to be completely
# out-grown by the other one.)

# ## Note: interfaces of constraint systems
#
# Internally, the community is connected via *interfaces*, which are small
# constraint trees (typically with no bounds attached) that describe parts of
# the constraint system that can be easily attached to other parts.
#
# The best kind of interface to choose generally differs from model to model.
# COBREXA gives a few "default" choices that cover a good part of sensible
# metabolic modeling. For example, if the model contains SBO annotations, we
# can ask for the interface created using the annotated reactions:
flux_balance_constraints(ecoli, interface = :sbo).interface

# If there are no annotations, we can still at least detect the boundary
# reactions and make an interface out of them:
flux_balance_constraints(ecoli, interface = :boundary).interface

# The default kind of interface in [`community_flux_balance_analysis`](@ref) is
# `:identifier_prefixes`, which relies on usual prefixes of reaction names
# (such as `EX_` for exchanges).

# Even if all of these methods fail, a suitable interface yourself can be
# produced manually. (Additionally, we can do useful stuff, such as removing
# the unnecessary bounds from the exchange descriptions.)
custom_model = flux_balance_constraints(ecoli)
custom_model *= remove_bounds(
    :interface^C.ConstraintTree(
        :biomass => custom_model.fluxes.BIOMASS_Ecoli_core_w_GAM,
        :exchanges => C.ConstraintTree(
            k => v for (k, v) in custom_model.fluxes if startswith(string(k), "EX_")
        ),
    ),
)
custom_model.interface.exchanges

# ## Connecting the community constraints manually
#
# To connect such interfaces into a community model, simply use function
# [`interface_constraints`](@ref) (which is how
# [`community_flux_balance_analysis`](@ref) constructs the community model
# internally via [`community_flux_balance_constraints`](@ref)). The assembly
# might look roughly as follows:

custom_community = interface_constraints(
    "bug1" => (
        custom_model * :handicap^C.Constraint(custom_model.fluxes.CYTBD.value, 0),
        0.2,
    ),
    "bug2" =>
        (custom_model * :handicap^C.Constraint(custom_model.fluxes.FBA.value, 0), 0.8),
    bound = r -> r == (:exchanges, :EX_glc__D_e) ? C.Between(-10, 0) : nothing,
)

custom_community.interface.exchanges

# For the model to work properly, we would need to add several other things,
# mainly the equal growth constraints (possibly via
# [`all_equal_constraints`](@ref)).
# [`community_flux_balance_constraints`](@ref) add these automatically, so we
# can equivalently just supply the constraint trees, and re-use the rest of the
# implementation:

custom_community = community_flux_balance_constraints(
    [
        "bug1" => (
            custom_model * :handicap^C.Constraint(custom_model.fluxes.CYTBD.value, 0),
            0.2,
        ),
        "bug2" => (
            custom_model * :handicap^C.Constraint(custom_model.fluxes.FBA.value, 0),
            0.8,
        ),
    ],
    ["EX_glc__D_e" => (-10.0, 0.0)],
)

# This can be solved with the usual means, reaching the same result as above:

custom_solution = optimized_values(
    custom_community,
    objective = custom_community.community_biomass.value,
    output = custom_community.community_biomass,
    optimizer = HiGHS.Optimizer,
)

@test isapprox(custom_solution, solution.community_biomass, atol = TEST_TOLERANCE) #src
