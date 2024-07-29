
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
# [`screen`](@ref) is a function that you can use to run many different
# simulations on many different variants of the model efficiently in parallel.

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import JSONFBCModels, HiGHS
import ConstraintTrees as C

model = load_model("e_coli_core.json")

# Screening is done as follows: We specify a range of values to examine (in
# the example below to an actual range of `-10` to `0`), but any vector-like
# list of things can be used), and write a short "analysis" function that takes
# one of the values from the range as an argument and runs an analysis for that
# given value.
#
# Here, we solve 11 different optimization problems for different bounds of the
# oxygen exchange:

screen(-10:0) do o2_bound
    c = flux_balance_constraints(model)
    c.fluxes.EX_o2_e.bound = C.EqualTo(o2_bound)
    optimized_values(
        c,
        objective = c.objective.value,
        optimizer = HiGHS.Optimizer,
        output = c.objective,
    )
end

# ## Screening in multiple dimensions

# The simplest way to screen through multi-dimensional arrays is to use an
# iterator product. In the result, we receive a whole matrix of results instead
# of a vector.

screen(Iterators.product(-10:2:0, 0:2:10)) do (o2_bound, co2_bound)
    c = flux_balance_constraints(model)
    c.fluxes.EX_o2_e.bound = C.EqualTo(o2_bound)
    c.fluxes.EX_co2_e.bound = C.EqualTo(co2_bound)
    optimized_values(
        c,
        objective = c.objective.value,
        optimizer = HiGHS.Optimizer,
        output = c.objective,
    )
end

# ## Screening through non-numeric properties

# If the problem at hand can not be expressed with mere ranges, we can specify
# vectors with any values. The following code examines the inter-dependency of
# blocking selected transport reactions together with selected exchanges, with
# 2 different bounds that implement the block. Because of 3-dimensional input,
# the result is expectably a 3-dimensional array:

screen(
    Iterators.product(
        [:H2Ot, :CO2t, :O2t, :NH4t], # a set of transport reactions
        [:EX_h2o_e, :EX_co2_e, :EX_o2_e, :EX_nh4_e], # a set of exchanges
        [C.Between(-0.1, 0.1), C.Between(-1, 1)], # bounds
    ),
) do (transport, exchange, bound)
    c = flux_balance_constraints(model)
    c.fluxes[transport].bound = 5 * bound
    c.fluxes[exchange].bound = 3 * bound
    optimized_values(
        c,
        objective = c.objective.value,
        optimizer = HiGHS.Optimizer,
        output = c.objective,
    )
end

# ## Annotating the results

# For reasons of tidyness, it is adviseable to annotate all values with their
# actual meaning directly in the arrays.
#
# With `screen`, the simplest (and usually sufficiently effective) way to do
# that is to return `Pair`s with annotation keys instead of plain values:

screen(
    Iterators.product(
        [:H2Ot, :CO2t, :O2t, :NH4t],
        [:EX_h2o_e, :EX_co2_e, :EX_o2_e, :EX_nh4_e],
        [C.Between(-0.1, 0.1), C.Between(-1, 1)],
    ),
) do (transport, exchange, bound)
    c = flux_balance_constraints(model)
    c.fluxes[transport].bound = 5 * bound
    c.fluxes[exchange].bound = 3 * bound
    (transport, exchange, bound.upper) => optimized_values(
        c,
        objective = c.objective.value,
        optimizer = HiGHS.Optimizer,
        output = c.objective,
    )
end

# Notably, this approach makes various indexing and ordering errors quite improbable.
