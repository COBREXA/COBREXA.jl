
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

Perform a Flux Variability Analysis (FVA) on the `model`, and return a
dictionary of flux ranges where the model is able to perform optimally.

The constraint system is constructed using [`flux_balance_constraints`](@ref),
and the variability is examined on all reaction's fluxes, or on the subset
given optionally in `reaction_subset` (e.g., `reaction_subset = ["PFK",
"ACALD"]`). The optimality tolerance can be specified with objective_bound
using e.g. [`relative_tolerance_bound`](@ref) or
[`absolute_tolerance_bound`](@ref); the default is 99% relative tolerance.

Parameter `workers` may be used to enable parallel or distributed processing;
the execution defaults to all available workers. Other parameters (esp.
`optimizer`) are internally forwarded to [`optimized_values`](@ref).

Use [`constraints_variability`](@ref) to customize the FVA execution.
"""
function flux_variability_analysis(
    model::A.AbstractFBCModel;
    objective_bound = relative_tolerance_bound(0.99),
    reaction_subset = nothing,
    optimizer,
    settings,
    workers = D.workers(),
)
    constraints = flux_balance_constraints(model)

    objective = constraints.objective_value

    objective_flux = optimized_values(
        constraints;
        objective = constraints.objective.value,
        output = constraints.objective,
        optimizer,
        settings,
    )

    isnothing(objective_flux) && return nothing

    constraints_variability(
        constraints *
        :objective_bound^C.Constraint(objective, objective_bound(objective_flux)),
        isnothing(reaction_subset) ? constraints.fluxes :
        let s = Set(Symbol.(reaction_subset))
            C.ConstraintTree(k => v for (k, v) in constraints.fluxes if k in s)
        end;
        optimizer,
        settings,
        workers,
    )
end

export flux_variability_analysis
