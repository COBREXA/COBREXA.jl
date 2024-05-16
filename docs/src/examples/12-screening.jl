
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

# # Screening through many model variants
#
# [`screen`](@ref) is a function that runs many model/simulation variants
# (ideally on an HPC) efficiently.

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import JSONFBCModels, GLPK
import ConstraintTrees as C

model = load_model("e_coli_core.json")

screen(-10:0) do o2_bound
    c = flux_balance_constraints(model)
    c.fluxes.EX_o2_e.bound = C.EqualTo(o2_bound)
    optimized_values(
        c,
        objective = c.objective.value,
        optimizer = GLPK.Optimizer,
        output = c.objective,
    )
end

# ## Screening in multiple dimensions

screen(Iterators.product(-10:2:0, 0:2:10)) do (o2_bound, co2_bound)
    c = flux_balance_constraints(model)
    c.fluxes.EX_o2_e.bound = C.EqualTo(o2_bound)
    c.fluxes.EX_co2_e.bound = C.EqualTo(co2_bound)
    optimized_values(
        c,
        objective = c.objective.value,
        optimizer = GLPK.Optimizer,
        output = c.objective,
    )
end

# ## Screening through non-numeric properties

screen(
    Iterators.product(
        [:H2Ot, :CO2t, :O2t, :NH4t], # a set of transport reactions
        [:EX_h2o_e, :EX_co2_e, :EX_o2_e, :EX_nh4_e], # a set of exchanges
        [C.Between(-0.1, 0.1), C.Between(-1, 1)],
    ),
) do (transport, exchange, bound)
    c = flux_balance_constraints(model)
    c.fluxes[transport].bound = 5 * bound
    c.fluxes[exchange].bound = 3 * bound
    optimized_values(
        c,
        objective = c.objective.value,
        optimizer = GLPK.Optimizer,
        output = c.objective,
    )
end

# ## Annotating the results

screen(
    Iterators.product(
        [:H2Ot, :CO2t, :O2t, :NH4t], # a set of transport reactions
        [:EX_h2o_e, :EX_co2_e, :EX_o2_e, :EX_nh4_e], # a set of exchanges
        [C.Between(-0.1, 0.1), C.Between(-1, 1)],
    ),
) do (transport, exchange, bound)
    c = flux_balance_constraints(model)
    c.fluxes[transport].bound = 5 * bound
    c.fluxes[exchange].bound = 3 * bound
    (transport, exchange, bound.upper) => optimized_values(
        c,
        objective = c.objective.value,
        optimizer = GLPK.Optimizer,
        output = c.objective,
    )
end
