# Copyright (c) 2025, University of Luxembourg                              #src
# Copyright (c) 2025, Heinrich-Heine University Duesseldorf                 #src
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

# # Steady community composition models

using COBREXA

# In the [community construction example](04-community-models.md), we have
# constructed a 2-member community of interacting *E. coli* with fixed
# abundances. Here we show how to explore the community abundances in a more
# dynamic way using [`community_composition_balance_constraints`](@ref), which
# follows the methodology established by the SteadyCom method.
#
# In short, we start with a similar community of knockouts, and determine what
# are the feasible ranges for the community abundances.
#
# As usual, we start with the toy *E. coli* model and a few packages:

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import JSONFBCModels
import HiGHS
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C

# We remove the artificial limit of glucose intake (so that we can move it to
# the community level later):
ecoli = load_model("e_coli_core.json", CM.Model)
ecoli.reactions["EX_glc__D_e"].lower_bound = -1000.0
ecoli.reactions["EX_glc__D_e"].upper_bound = 1000.0

# We create a few reaction knockouts:
knockouts = ["CYTBD", "FBA"];
ko_ecolis = Dict(begin
    m = deepcopy(ecoli)
    m.reactions[r].lower_bound = m.reactions[r].upper_bound = 0
    Symbol(r) => m
end for r in knockouts)

# ## Finding feasible community compositions
#
# Now, let's use [`community_composition_variability_analysis`](@ref) to
# determine the feasible range of abundances if the community if these 2
# reaction knockouts is forced to grow 0.5 g/gDW/h:

res = community_composition_variability_analysis(
    ko_ecolis,
    0.5,
    ["EX_glc__D_e" => (-10.0, 0.0)],
    optimizer = HiGHS.Optimizer,
)

# It seems like the FBA-knockout has to be present for the community to be
# feasible at this rate.

@test isapprox(first(res.CYTBD), 0, atol = TEST_TOLERANCE) #src
@test isapprox(last(res.CYTBD), 0.6642931585890363, atol = TEST_TOLERANCE) #src
@test isapprox(first(res.FBA), 0.33570684141096374, atol = TEST_TOLERANCE) #src
@test isapprox(last(res.FBA), 1, atol = TEST_TOLERANCE) #src

# ## Finding the community composition at optimum growth
#
# The variable community composition allows us to also scan for optimal growth
# (via simple interval splitting), determining the community composition near
# this optimum:
res = community_composition_balance_analysis(
    ko_ecolis,
    1.0,
    ["EX_glc__D_e" => (-10.0, 0.0)],
    optimizer = HiGHS.Optimizer,
    tolerance = 1e-4,
)

# The achieved growth is reported in the result tree:
res.growth

@test isapprox(res.growth, 0.7039794921875, atol = TEST_TOLERANCE) #src

# We can also have a look at the abundances at optimum. In this case, it looks
# like CYTBD knockout will be practically missing (within the tolerance
# requested above, it is technically zero):
res.abundances

@test isapprox(sum(values(res.abundances)), 1.0, atol = TEST_TOLERANCE) #src
@test isapprox(res.abundances.FBA, 1.0, atol = 1e-4) #src
