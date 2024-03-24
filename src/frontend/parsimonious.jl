
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

A constraint system like from [`flux_balance_constraints`](@ref), but with the
parsimonious objective present under key `parsimonious_objective`. Best used
via [`parsimonious_flux_balance_analysis`](@ref).
"""
function parsimonious_flux_balance_constraints(model::A.AbstractFBCModel)
    constraints = flux_balance_constraints(model)

    constraints *
    :parsimonious_objective^C.Constraint(squared_sum_value(constraints.fluxes))
end

export parsimonious_flux_balance_constraints

"""
$(TYPEDSIGNATURES)

Compute a parsimonious flux solution for the `model`, using the constraints
given by [`parsimonious_flux_balance_constraints`](@ref).

In short, the objective value of the parsimonious solution should be the same
as the one from [`flux_balance_analysis`](@ref), except the squared sum of
reaction fluxes is minimized. If there are multiple possible fluxes that
achieve a given objective value, parsimonious flux thus represents the "minimum
energy" one, which is arguably more realistic.

Solver configuration arguments are forwarded to
[`parsimonious_optimized_values`](@ref).
"""
parsimonious_flux_balance_analysis(
    model::A.AbstractFBCModel;
    tolerances = relative_tolerance_bound.(1 .- [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]),
    kwargs...,
) = frontend_parsimonious_optimized_values(
    parsimonious_flux_balance_constraints,
    model;
    objective = x -> x.objective.value,
    parsimonious_objective = x -> x.parsimonious_objective.value,
    tolerances,
    kwargs...,
)

export parsimonious_flux_balance_analysis

"""
$(TYPEDSIGNATURES)

Like [`parsimonious_flux_balance_constraints`](@ref), but uses a L1 metric for
solving the parsimonious problem. The `parsimonious_objective` constraint is
thus linear.
"""
function linear_parsimonious_flux_balance_constraints(model::A.AbstractFBCModel; kwargs...)
    constraints = flux_balance_constraints(model)
    constraints =
        constraints +
        :fluxes_forward^unsigned_positive_contribution_variables(constraints.fluxes) +
        :fluxes_reverse^unsigned_negative_contribution_variables(constraints.fluxes)
    return constraints *
           :directional_flux_balance^sign_split_constraints(
               positive = constraints.fluxes_forward,
               negative = constraints.fluxes_reverse,
               signed = constraints.fluxes,
           ) *
           :parsimonious_objective^C.Constraint(
               sum_value(constraints.fluxes_forward, constraints.fluxes_reverse),
           )
end

export linear_parsimonious_flux_balance_constraints

"""
$(TYPEDSIGNATURES)

Like [`parsimonious_flux_balance_analysis`](@ref), but uses the L1-metric
parsimonious system given by
[`linear_parsimonious_flux_balance_constraints`](@ref).

In turn, the solution is often faster, does not require a solver capable of
quadratic objectives, and has many beneficial properties of the usual
parsimonious solutions (such as the general lack of unnecessary loops). On the
other hand, like with plain flux balance analysis there is no strong guarantee
of uniqueness of the solution.

Solver configuration arguments are forwarded to
[`parsimonious_optimized_values`](@ref).
"""
linear_parsimonious_flux_balance_analysis(
    model::A.AbstractFBCModel;
    tolerances = relative_tolerance_bound.(1 .- [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]),
    kwargs...,
) = frontend_parsimonious_optimized_values(
    linear_parsimonious_flux_balance_constraints,
    model;
    objective = x -> x.objective.value,
    parsimonious_objective = x -> x.parsimonious_objective.value,
    tolerances,
    kwargs...,
)

export linear_parsimonious_flux_balance_analysis
