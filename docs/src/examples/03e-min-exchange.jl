
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

# # Parsimonious growth medium
#
# [Parsimonious flux balance analysis](03b-parsimonious-flux-balance.md) may be
# easily combined with the [unidirectional reaction splitting
# tools](03d-unidirectional.md) to create an analysis that minimizes the total
# nutrient uptake flux. This is beneficial as a complementary, less "refined"
# but much more computationally efficient alternative to the [growth medium
# optimization](05h-medium.md) (which requires a demanding mixed-integer
# programming).
#
# For the demonstration, we use the toy E. coli model:

using COBREXA
import ConstraintTrees as C #src

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import JSONFBCModels, HiGHS
model = load_model("e_coli_core.json")

# First, we create a model with extra direction-splitted variables and
# constraints for the negative exchange reactions:

c = flux_balance_constraints(model, interface = :identifier_prefixes)
c += :in_exchanges^unsigned_negative_contribution_variables(c.interface.exchanges)
c *=
    :out_exchanges^unsigned_positive_contribution_constraints(
        c.interface.exchanges,
        c.in_exchanges,
    )

# For the parsimonious optimization, we need a good guess on the objective
# value that we will require; which we obtain from an optimal solution of the
# original model:

objective_value = flux_balance_analysis(model, optimizer = HiGHS.Optimizer).objective

# Parsimonious analyses need a clear specification of a solution tolerance; we
# use a simple absolute difference bound:

tolerances = [absolute_tolerance_bound(0.001)]

# We can now run a linear "parsimonious intake" analysis as follows:
result = parsimonious_optimized_values(
    c;
    objective = c.objective.value,
    parsimonious_objective = sum_value(c.in_exchanges),
    objective_value,
    tolerances,
    optimizer = HiGHS.Optimizer,
    output = c.in_exchanges,
)

# If suitable, one may also choose to run a quadratic (L2) parsimonious
# analysis, simply by swapping the objective:
l2_result = parsimonious_optimized_values(
    c;
    objective = c.objective.value,
    parsimonious_objective = sum_value(c.in_exchanges),
    objective_value,
    tolerances,
    optimizer = HiGHS.Optimizer,
    output = c.in_exchanges,
)

@test all(values(C.zip(result, l2_result, Bool) do a, b #src
    isapprox(a, b, atol = TEST_TOLERANCE) #src
end)) #src

# Because of the model simplicity, the solution coincides with the linear one.
#
# To better demonstrate the difference between L1 and L2 solutions that more
# complex models would exhibit, we arbitrarily choose to relax the tolerance
# (indirectly allowing non-unique solution):

tolerances = [relative_tolerance_bound(0.5)]

# The linear solution is appropriately reduced:
parsimonious_optimized_values(
    c;
    objective = c.objective.value,
    parsimonious_objective = sum_value(c.in_exchanges),
    objective_value,
    tolerances,
    optimizer = HiGHS.Optimizer,
    output = c.in_exchanges,
)

# The quadratic solution additionally reduces the "extremes" in the flux
# profile, since the quadratic cost makes them unfavorable:
parsimonious_optimized_values(
    c;
    objective = c.objective.value,
    parsimonious_objective = squared_sum_value(c.in_exchanges),
    objective_value,
    tolerances,
    optimizer = HiGHS.Optimizer,
    output = c.in_exchanges,
)
