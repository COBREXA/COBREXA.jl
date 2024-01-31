
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

Find a solution of the "minimization of metabolic adjustment" (MOMA) analysis
for the `model`, which is the "closest" feasible solution to the given
`reference_fluxes`, in the sense of squared-sum distance. The minimized
squared distance (the objective) is present in the result tree as
`minimal_adjustment_objective`.

This is often used for models with smaller feasible region than the reference
models (typically handicapped by a knockout, nutritional deficiency or a
similar perturbation). MOMA solution then gives an expectable "easiest"
adjustment of the organism towards a somewhat working state.

Reference fluxes that do not exist in the model are ignored (internally, the
objective is constructed via [`squared_sum_error_value`](@ref)).

Additional parameters are forwarded to [`optimized_values`](@ref).
"""
function metabolic_adjustment_minimization_constraints(
    model::A.AbstractFBCModel,
    reference_fluxes::Dict{Symbol,Float64},
)
    constraints = flux_balance_constraints(model)
    constraints *
    :minimal_adjustment_objective^C.Constraint(
        squared_sum_error_value(constraints.fluxes, x -> get(reference_fluxes, x, nothing)),
    )
end

"""
$(TYPEDSIGNATURES)

A slightly easier-to-use version of
[`metabolic_adjustment_minimization_constraints`](@ref) that computes the
reference flux as the parsimonious optimal solution of the `reference_model`.
The reference flux is calculated using `reference_optimizer` and
`reference_modifications`, which default to the `optimizer` and `settings`.

Other arguments are forwarded to the internal call of
[`parsimonious_optimized_values`](@ref).

Returns `nothing` if no feasible solution is found.
"""
function metabolic_adjustment_minimization_constraints(
    model::A.AbstractFBCModel,
    reference_model::A.AbstractFBCModel;
    kwargs...,
)
    reference_constraints = parsimonious_flux_balance_constraints(reference_model)
    reference_fluxes = parsimonious_optimized_values(
        reference_constraints;
        objective = reference_constraints.objective.value,
        parsimonious_objective = reference_constraints.parsimonious_objective.value,
        output = reference_constraints.fluxes,
        kwargs...,
    )
    isnothing(reference_fluxes) && return nothing
    metabolic_adjustment_minimization_constraints(model, reference_fluxes)
end

export metabolic_adjustment_minimization_constraints

"""
$(TYPEDSIGNATURES)

TODO
"""
metabolic_adjustment_minimization_analysis(model::A.AbstractFBCModel, args...; kwargs...) =
    frontend_optimized_values(
        metabolic_adjustment_minimization_constraints,
        model,
        args...;
        objective = x -> x.minimal_adjustment_objective.value,
        sense = Minimal,
        kwargs...,
    )

export metabolic_adjustment_minimization_analysis
#
# TODO this needs to pass in the reference_optimizer and other args properly (there are actually 3 sets of optimizer args: reference 1st problem, reference parsimonious problem, and MOMA problem.

"""
$(TYPEDSIGNATURES)

Like [`metabolic_adjustment_minimization_constraints`](@ref) but optimizes the L1 norm.
This typically produces a sufficiently good result with less resources,
depending on the situation. See documentation of
[`linear_parsimonious_flux_balance_analysis`](@ref) for some of the
considerations.
"""
function linear_metabolic_adjustment_minimization_constraints(
    model::A.AbstractFBCModel,
    reference_fluxes::Dict{Symbol,Float64},
)
    constraints = flux_balance_constraints(model)

    difference = C.zip(ct.fluxes, C.Tree(reference_fluxes)) do orig, ref
        C.Constraint(orig.value - ref)
    end

    difference_split_variables =
        C.variables(keys = keys(difference), bounds = C.Between(0, Inf))
    constraints += :reference_positive_diff^difference_split_variables
    constraints += :reference_negative_diff^difference_split_variables

    # `difference` actually doesn't need to go to the CT, but we include it
    # anyway to calm the curiosity of good neighbors.
    constraints *
    :reference_diff^difference *
    :reference_directional_diff_balance^sign_split_constraints(
        positive = constraints.reference_positive_diff,
        negative = constraints.reference_negative_diff,
        signed = difference,
    ) *
    :linear_minimal_adjustment_objective^C.Constraint(
        sum_value(constraints.reference_positive_diff, constraints.reference_negative_diff),
    )
end

"""
$(TYPEDSIGNATURES)

Like [`metabolic_adjustment_minimization_constraints`](@ref) but optimizes the L1 norm.
This typically produces a sufficiently good result with less resources,
depending on the situation. See documentation of
[`linear_parsimonious_flux_balance_analysis`](@ref) for some of the
considerations.
"""
function linear_metabolic_adjustment_minimization_constraints(
    model::A.AbstractFBCModel,
    reference_model::A.AbstractFBCModel;
    kwargs...,
)
    reference_constraints = parsimonious_flux_balance_constraints(reference_model)
    reference_fluxes = parsimonious_optimized_values(
        reference_constraints;
        objective = reference_constraints.objective.value,
        parsimonious_objective = reference_constraints.parsimonious_objective.value,
        output = reference_constraints.fluxes,
        kwargs...,
    )
    isnothing(reference_fluxes) && return nothing

    linear_metabolic_adjustment_minimization_constraints(model, reference_fluxes)
end

export linear_metabolic_adjustment_minimization_constraints

"""
$(TYPEDSIGNATURES)

TODO
"""
linear_metabolic_adjustment_minimization_analysis(
    model::A.AbstractFBCModel,
    args...;
    kwargs...,
) = frontend_optimized_values(
    linear_metabolic_adjustment_minimization_constraints,
    model,
    args...;
    objective = x -> x.linear_minimal_adjustment_objective.value,
    sense = Minimal,
    kwargs...,
)

export linear_metabolic_adjustment_minimization_analysis

#TODO: somehow forward optimizer args to linear moma constraints builder (see above)
