
# Copyright (c) 2025, University of Luxembourg
# Copyright (c) 2025, Heinrich-Heine University Duesseldorf
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

Make a medium-optimizing mixed-integer system from a FBC model.

The medium is assumed to consist of intakes (negative fluxes) in the selected
`exchange_reactions`. If these are not specified, the model must contain
exchanges that are identifiable as an exchange interface by
[`flux_balance_constraints`](@ref) to be used as the medium (by default, these
are identified from identifier prefixes; parameters `interface` and
`interface_name` are forwarded).

The output system will be constrainted to reach at least `objective_target`
flux through the objective function.

`exchange_reaction_cost` should assign a numeric cost of inclusion of each of
the reactions in the `universal_model`; by default all are assigned equal cost
of `1`. `max_cost` puts an optional maximum limit on the cost, which may help
the solver to avoid exploring unnecessarily complex solutions.  `known_flags`
may supply previous solutions of the same system; these will be made infeasible
in the output constraint system in order to allow discovery of different ones.

Additional arguments are forwarded to `flux_balance_constraints` that converts
`model` to constraints.
"""
function medium_optimization_constraints(
    model::A.AbstractFBCModel,
    objective_target::Float64;
    exchange_reactions::Maybe{Vector{String}} = nothing,
    exchange_reaction_cost = _ -> 1.0,
    max_cost = Inf,
    known_flags::Vector{C.Tree{Float64}} = C.Tree{Float64}[],
    interface::Maybe{Symbol} = :identifier_prefixes,
    interface_name = :interface,
    kwargs...,
)
    m = flux_balance_constraints(model; interface, interface_name, kwargs...)
    m.objective.bound = C.Between(objective_target, Inf)

    medium_optimization_constraints(
        m,
        isnothing(exchange_reactions) ? m[interface_name].exchanges :
        C.ConstraintTree(
            Symbol(str) => m.fluxes[Symbol(str)] for str in exchange_reactions
        );
        exchange_cost = sym -> exchange_reaction_cost(String(sym)),
        max_cost,
        known_flags,
    )
end

"""
$(TYPEDSIGNATURES)

Make a medium optimization system out of `constraints` that contain
`exchanges`. In the system, the total cost of exchanges (optionally weighted by
`exchange_cost` function) is minimized and limited by `max_cost`.

`known_flags` may contain a vector of `medium_flags` subtrees from the previous
solutions; solutions with the same flags as in the vector are avoided.
"""
function medium_optimization_constraints(
    constraints::C.ConstraintTree,
    exchanges::C.ConstraintTree;
    exchange_cost = _ -> 1.0,
    max_cost = Inf,
    known_flags::Vector{C.Tree{Float64}} = C.Tree{Float64}[],
)
    joined =
        C.ConstraintTree(:system => constraints, :exchanges => exchanges) +
        :medium_flags^C.variables_for(exchanges) do _
            Switch(0, 1)
        end

    return C.ConstraintTree(
        :system => joined.system,
        :exchange_flag_bounds => C.zip(joined.exchanges, joined.medium_flags) do x, flag
            value_scaled_bound_constraint(
                x.value,
                greater_or_equal_interval(x.bound),
                flag.value,
            )
        end,
        :medium_flags => joined.medium_flags,
        (
            Symbol(:known_medium_, i) =>
                gap_filling_known_fill_constraint(joined.medium_flags, kf) for
            (i, kf) in enumerate(known_flags)
        )...,
        :medium_size => C.Constraint(
            C.sum((v.value for (k, v) in joined.medium_flags), init = zero(C.LinearValue)),
        ),
        :medium_cost => C.Constraint(
            C.sum(
                (exchange_cost(k) * v.value for (k, v) in joined.medium_flags),
                init = zero(C.LinearValue),
            ),
            C.Between(0, max_cost),
        ),
    )
end

export medium_optimization_constraints

"""
$(TYPEDSIGNATURES)

Run a medium optimization analysis on a constraint system as built by
[`medium_optimization_constraints`](@ref).
"""
medium_optimization_analysis(args...; kwargs...) = frontend_optimized_values(
    medium_optimization_constraints,
    args...;
    objective = x -> x.medium_cost.value,
    sense = Minimal,
    kwargs...,
)

export medium_optimization_analysis
