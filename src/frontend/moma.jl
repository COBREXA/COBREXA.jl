
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

Keyword arguments are discarded for compatibility with the other overload.
"""
function metabolic_adjustment_minimization_constraints(
    model::A.AbstractFBCModel,
    reference_fluxes::C.Tree;
    _...,
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

Find a solution of the "minimization of metabolic adjustment" (MOMA) analysis
for the `model`, which is the "closest" feasible solution to the solution given
in the second argument, which is either `reference_fluxes` or `reference_model`
(see documentation of [`metabolic_adjustment_minimization_constraints`](@ref)),
in the sense of squared-sum distance. The minimized squared distance (the
objective) is present in the result tree as `minimal_adjustment_objective`.

If the second argument is a reference model, it is solved using a
[`parsimonious_flux_balance_analysis`](@ref) with the optimizer and settings
parameters for the 2 steps set by keyword arguments prefixed by `reference_`.

This is often used for models with smaller feasible region than the reference
models (typically handicapped by a knockout, a nutritional deficiency or a
similar perturbation). MOMA solution then gives an expectable "easiest"
adjustment of the organism towards a somewhat working state.

Reference fluxes that do not exist in the `model` are ignored (internally, the
objective is constructed via [`squared_sum_error_value`](@ref)).
"""
metabolic_adjustment_minimization_analysis(
    model::A.AbstractFBCModel,
    args...;
    optimizer,
    settings = [],
    reference_parsimonious_optimizer = optimizer,
    reference_parsimonious_settings = settings,
    reference_optimizer = optimizer,
    reference_settings = settings,
    kwargs...,
) = frontend_optimized_values(
    metabolic_adjustment_minimization_constraints,
    model,
    args...;
    builder_kwargs = (
        optimizer = reference_optimizer,
        settings = reference_settings,
        parsimonious_optimizer = reference_parsimonious_optimizer,
        parsimonious_settings = reference_parsimonious_settings,
    ),
    objective = x -> x.minimal_adjustment_objective.value,
    sense = Minimal,
    optimizer,
    settings,
    kwargs...,
)

export metabolic_adjustment_minimization_analysis

"""
$(TYPEDSIGNATURES)

Like [`metabolic_adjustment_minimization_constraints`](@ref) but optimizes the
L1 distance from `reference_fluxes`.

Keyword arguments are discarded for compatibility with the other overload.
"""
function linear_metabolic_adjustment_minimization_constraints(
    model::A.AbstractFBCModel,
    reference_fluxes::C.Tree;
    _...,
)
    constraints = flux_balance_constraints(model)

    # `difference` actually doesn't need to go to the CT, but we include it
    # anyway to calm the curiosity of good neighbors.
    constraints *= :reference_diff^C.zip(constraints.fluxes, reference_fluxes) do orig, ref
        C.Constraint(orig.value - ref)
    end

    # split the reference_diff into positive and negative contributions
    constraints +=
        :reference_positive_diff^unsigned_positive_contribution_variables(
            constraints.reference_diff,
        )
    constraints *=
        :reference_negative_diff^unsigned_negative_contribution_constraints(
            constraints.reference_diff,
            constraints.reference_positive_diff,
        )

    constraints *
    :linear_minimal_adjustment_objective^C.Constraint(
        sum_value(constraints.reference_positive_diff, constraints.reference_negative_diff),
    )
end

"""
$(TYPEDSIGNATURES)

Like [`metabolic_adjustment_minimization_constraints`](@ref) but the output
constraints optimize the L1 distance from the linear-parsimonious solution of
the `reference_model`.
"""
function linear_metabolic_adjustment_minimization_constraints(
    model::A.AbstractFBCModel,
    reference_model::A.AbstractFBCModel;
    kwargs...,
)
    reference_constraints = linear_parsimonious_flux_balance_constraints(reference_model)
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

Perform a linear minimization of metabolic adjustment analysis (l-MOMA) on
`model`. The reference is given by the second argument, which is either a
`reference_flux` or a `reference_model` (the second argument is forwarded to
[`linear_metabolic_adjustment_minimization_constraints`](@ref)).

While the solution is "less uniquely defined" than with fully quadratic
[`metabolic_adjustment_minimization_analysis`](@ref), the linear variant
typically produces a sufficiently good result with much less resources. See
documentation of [`linear_parsimonious_flux_balance_analysis`](@ref) for some
of the considerations.
"""
linear_metabolic_adjustment_minimization_analysis(
    model::A.AbstractFBCModel,
    args...;
    optimizer,
    settings = [],
    reference_parsimonious_optimizer = optimizer,
    reference_parsimonious_settings = settings,
    reference_optimizer = optimizer,
    reference_settings = settings,
    kwargs...,
) = frontend_optimized_values(
    linear_metabolic_adjustment_minimization_constraints,
    model,
    args...;
    builder_kwargs = (
        optimizer = reference_optimizer,
        settings = reference_settings,
        parsimonious_optimizer = reference_parsimonious_optimizer,
        parsimonious_settings = reference_parsimonious_settings,
    ),
    objective = x -> x.linear_minimal_adjustment_objective.value,
    sense = Minimal,
    optimizer,
    settings,
    kwargs...,
)

export linear_metabolic_adjustment_minimization_analysis
