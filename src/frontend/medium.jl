
# Copyright (c) 2024, University of Luxembourg
# Copyright (c) 2024, Heinrich-Heine University Duesseldorf
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
        :exchange_fill_bounds => C.zip(joined.exchanges, joined.medium_flags) do x, b
            value_scaled_bound_constraint(x.value, x.bound, b.value)
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

medium_optimization_analysis(args...; kwargs...) = frontend_optimized_values(
    medium_optimization_constraints,
    args...;
    objective = x -> x.medium_cost.value,
    sense = Minimal,
    kwargs...,
)

export medium_optimization_analysis
