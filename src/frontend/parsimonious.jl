
# Copyright (c) 2021-2024, University of Luxembourg
# Copyright (c) 2021-2024, Heinrich-Heine University Duesseldorf
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
$(TYPEDSIGNATURES)

Compute a parsimonious flux solution for the given `model`. In short, the
objective value of the parsimonious solution should be the same as the one from
[`flux_balance_analysis`](@ref), except the squared sum of reaction fluxes is minimized.
If there are multiple possible fluxes that achieve a given objective value,
parsimonious flux thus represents the "minimum energy" one, thus arguably more
realistic. The optimized squared distance is present in the result as
`parsimonious_objective`.

Most arguments are forwarded to [`parsimonious_optimized_values`](@ref),
with some (objectives) filled in automatically to fit the common processing of
FBC models, and some (`tolerances`) provided with more practical defaults.

Similarly to the [`flux_balance_analysis`](@ref), returns a tree with the optimization
solutions of the shape as given by [`flux_balance_constraints`](@ref).
"""
function parsimonious_flux_balance_analysis(
    model::A.AbstractFBCModel,
    optimizer;
    tolerances = relative_tolerance_bound.(1 .- [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]),
    kwargs...,
)
    constraints = flux_balance_constraints(model)
    parsimonious_objective = squared_sum_value(constraints.fluxes)
    parsimonious_optimized_values(
        constraints * :parsimonious_objective^C.Constraint(parsimonious_objective);
        optimizer,
        objective = constraints.objective.value,
        parsimonious_objective,
        tolerances,
        kwargs...,
    )
end

export parsimonious_flux_balance_analysis

"""
$(TYPEDSIGNATURES)

Like [`parsimonious_flux_balance_analysis`](@ref), but uses a L1 metric for
solving the parsimonious problem.

In turn, the solution is often faster, does not require a solver capable of
quadratic objectives, and has many beneficial properties of the usual
parsimonious solutions (such as the general lack of unnecessary loops). On the
other hand, like with plain flux balance analysis there is no strong guarantee
of uniqueness of the solution.
"""
function linear_parsimonious_flux_balance_analysis(
    model::A.AbstractFBCModel,
    optimizer;
    tolerances = relative_tolerance_bound.(1 .- [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]),
    kwargs...,
)
    constraints = flux_balance_constraints(model)
    constraints =
        constraints +
        :fluxes_forward^unsigned_positive_contribution_variables(ct.fluxes) +
        :fluxes_reverse^unsigned_negative_contribution_variables(ct.fluxes)
    constraints *=
        :directional_flux_balance^sign_split_constraints(
            positive = ct.fluxes_forward,
            negative = ct.fluxes_reverse,
            signed = ct.fluxes,
        )

    parsimonious_objective = sum_value(ct.fluxes_forward, ct.fluxes_reverse)

    parsimonious_optimized_values(
        constraints * :parsimonious_objective^C.Constraint(parsimonious_objective);
        optimizer,
        objective = constraints.objective.value,
        parsimonious_objective,
        tolerances,
        kwargs...,
    )
end

export linear_parsimonious_flux_balance_analysis