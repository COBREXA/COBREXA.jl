
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

# # Growth medium optimization
#
# With mixed-integer linear programming, it is possible to find the "simplest"
# substrate medium that can support the growth of the model. This
# medium-optimization analysis minimizes the total amount of exchange reactions
# that are active (i.e., the flux is not zero). Medium optimization can be used
# to scan for diverse feeds that can support your model, and to identify
# metabolites that are crucial for growth.
#
# Here, we will use the toy *E. coli* toy model and see what are the diverse
# "minimal" nutrient combinations that can support its growth.
#
#md # !!! warning "Minimal vs. smallest"
#md #     Medium optimization identifies a set of exchange reaction of minimal *size*, i.e., the smallest possible count of exchange reactions active. This does not correspond to the smallest or most efficient exchange flux that is achievable -- to identify minimal fluxes, use [parsimonious analysis](03b-parsimonious-flux-balance.md).

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import JSONFBCModels, HiGHS
model = load_model("e_coli_core.json")

# To get a good guess on how much the model can grow on an arbitrary medium, we
# solve it using FBA:

fba_result = flux_balance_analysis(model; optimizer = HiGHS.Optimizer);
fba_result.objective

# ## Find a minimal set of exchanges

# Medium optimization is implemented in [`medium_optimization_analysis`](@ref).
#
# The function needs to know a way to identify the exchange reactions; by
# default it uses [`flux_balance_constraints`](@ref) with parameter `interface
# = :identifier_prefixes`, and finds the exchanges via the generated interface.
# If one aims to optimize a different set of exchanges, it is possible to
# specify these via argument `exchange_reactions`.

x = medium_optimization_analysis(
    model,
    fba_result.objective * 0.999,
    optimizer = HiGHS.Optimizer,
)

# We can find the medium from the "flag" variables:

medium = [k for (k, v) in x.medium_flags if v != 0]

@test length(medium) == 7 #src

y = medium_optimization_analysis(
    model,
    fba_result.objective * 0.1,
    optimizer = HiGHS.Optimizer,
    known_flags = [x.medium_flags],
)

# TODO limit "medium" to intakes only
