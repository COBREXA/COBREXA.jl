
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

# # Split unidirectional reactions in models
#
# By default, the constraint system produced by e.g.
# [`flux_balance_constraints`](@ref) assumes a single variable for each
# reaction flux that may be both positive and negative (depending on the
# reaction). This example explains several ways to "split" such bidirectional
# fluxes into unidirectional "forward" and "reverse" parts. This is useful in
# modeling of capacity constraints (such system can be found e.g. in
# [`enzyme_constrained_flux_balance_analysis`](@ref) and
# [`linear_parsimonious_flux_balance_analysis`](@ref)) and many other
# purposes.
#
# Here we show how to create such system for the toy E. coli model:

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import JSONFBCModels
import HiGHS
import ConstraintTrees as C

model = load_model("e_coli_core.json")

# We will only work with the constraint representation of the model. The fluxes
# there are bidirectional:
c = flux_balance_constraints(model);
c.fluxes

# As the simplest approach, we can allocate 2 sets of variables for the forward
# and reverse fluxes via [`sign_split_variables`](@ref) and connect them to the
# fluxes with [`sign_split_constraints`](@ref). These functions ensure that the
# bounds of the unidirectional fluxes are within the expectable limit, and,
# respectively, that the original fluxes are constrained to match the sum of
# the split directions:

c += sign_split_variables(c.fluxes, positive = :fluxes_forward, negative = :fluxes_reverse);
c *=
    :directional_flux_balance^sign_split_constraints(
        positive = c.fluxes_forward,
        negative = c.fluxes_reverse,
        signed = c.fluxes,
    )

# We can solve this system as usual and observe the results as usual

x = optimized_values(c, objective = c.objective.value, optimizer = HiGHS.Optimizer)

C.zip(tuple, x.fluxes, x.fluxes_forward, x.fluxes_reverse, Tuple)

@test all( #src
    isapprox(0, atol = TEST_TOLERANCE), #src
    values( #src
        C.zip( #src
            (f, p, n) -> (p - n - f), #src
            x.fluxes, #src
            x.fluxes_forward, #src
            x.fluxes_reverse, #src
            Float64, #src
        ), #src
    ), #src
) #src

# ## Simplifying the system using asymmetric construction
#
# If used naively, the above construction uses 3 variables and many constraints
# for each flux, which is not quite efficient. To ameliorate the usage of
# solver resources, one may construct a slightly simpler (but asymmetric)
# system that only uses 2 variables:

c2 = flux_balance_constraints(model);
c2 += :fluxes_forward^unsigned_positive_contribution_variables(c2.fluxes);
c2 *=
    :fluxes_reverse^unsigned_negative_contribution_constraints(c2.fluxes, c2.fluxes_forward);

# This way, only forward fluxes are allocated as variables, and reverse fluxes
# are "computed" as linearly dependent values. Additionally, since the bounds
# on the forward and reverse fluxes completely subsume the original bounds on
# fluxes, one can further simplify the system by removing the original bounds:

c2.fluxes = remove_bounds(c2.fluxes)

# The system solves just like the "symmetric" one:
x2 = optimized_values(c2, objective = c2.objective.value, optimizer = HiGHS.Optimizer)

@test isapprox(x.objective, x2.objective, atol = TEST_TOLERANCE) #src

# We can also compare the raw variable counts:

(C.var_count(c), C.var_count(c2))

@test C.var_count(c) == 285 #src
@test C.var_count(c2) == 190 #src

# ## Simplifying the system by removing original variables
#
# If one can assume that the original system is just allocated variables with
# no other semantics attached, one can reduce the variable and constraint count
# even in the "nicer" symmetric case from above.
#
# In particular, it is possible to substitute a combination of forward and
# reverse flux for the bidirectional variables, which allows them to be pruned
# out of the system together with their original associated bounds:

subst_vals = [C.variable(; idx).value for idx = 1:C.var_count(c)]

c.fluxes = C.zip(c.fluxes, c.fluxes_forward, c.fluxes_reverse) do f, p, n
    (var_idx,) = f.value.idxs
    subst_value = p.value - n.value
    subst_vals[var_idx] = subst_value
    C.Constraint(subst_value) # the bidirectional bound is dropped here
end

c = C.prune_variables(C.substitute(c, subst_vals))

# The variable count is now reduced, and the system again solves just as the
# original:

C.var_count(c)
@test C.var_count(c) == 190 #src

#

x = optimized_values(c, objective = c.objective.value, optimizer = HiGHS.Optimizer);
x.objective

@test isapprox(x.objective, x2.objective, atol = TEST_TOLERANCE) #src

# The bidirectional flux values are computed transparently in the result:

x.fluxes
