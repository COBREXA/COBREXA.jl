
# Copyright (c) 2021-2025, University of Luxembourg
# Copyright (c) 2021-2025, Heinrich-Heine University Duesseldorf
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

Construct an instance of a linear community Flux Balance Analysis (cFBA) model.
The relative abundances of the organisms are known in advance; this function
predicts the maximum achievable community growth rate.

`model_abundances` is a dictionary-like object that maps model identifiers to
tuples of models (usually subtypes of `AbstractFBCModel`) and their abundances
(such as: `"bug1" => (bug1, 0.5), ...`). `community_exchange_bounds` is a
dictionary-like object that can be additionally used to restrict selected
community exchange reactions (keys should be reaction IDs, the values are
converted to `ConstraintTrees`-like bounds). Bounds otherwise default to
parameter `default_community_exchange_bound`, which itself defaults to
`nothing` (i.e., unbounded).

If required, constraint trees may be supplied instead of `AbstracFBCModel`s in
`model_abundances`. These must provide an interface compatible with
`interface_exchanges` and `interface_biomass`.

`interface` is forwarded to [`flux_balance_constraints`](@ref).
`interface_exchanges` and `interface_biomass` are used to pick up the correct
interface part to contribute to the community exchanges and community biomass.
`wrap` is forwarded to internal call of [`interface_constraints`](@ref).
"""
function community_flux_balance_constraints(
    model_abundances,
    community_exchange_bounds = Dict();
    interface = :identifier_prefixes,
    interface_exchanges = x -> x.interface.exchanges,
    interface_biomass = x -> x.interface.biomass,
    default_community_exchange_bound = nothing,
    wrap = identity,
)
    @assert length(model_abundances) >= 1 "at least one community member is required"
    @assert isapprox(sum(a for (_, (_, a)) in model_abundances), 1) "community member abundances must sum to 1"

    bounds_lookup = Dict(community_exchange_bounds)

    make_member(k, m, _) = throw(DomainError(m, "unsupported community member type ($k)"))
    make_member(k, m::A.AbstractFBCModel, a) =
        make_member(k, flux_balance_constraints(m; interface), a)
    make_member(_, c::C.ConstraintTree, a) = (c, interface_exchanges(c), a)

    constraints = interface_constraints(
        (Symbol(k) => make_member(k, m, a) for (k, (m, a)) in model_abundances);
        out_interface = :community_exchanges,
        out_balance = :community_balance,
        bound = x ->
            get(bounds_lookup, String(last(x)), default_community_exchange_bound),
        wrap,
    )

    growth_sums = [
        Symbol(k) => C.Constraint(sum_value(interface_biomass(constraints[Symbol(k)])))
        for (k, _) in model_abundances
    ]

    # TODO: this should be prefixed by `:community` to avoid name collisions.
    constraints *
    :equal_growth^all_equal_constraints(
        last(growth_sums[1]).value,
        C.ConstraintTree(growth_sums[2:end]),
    ) *
    :community_biomass^last(growth_sums[1])
end

export community_flux_balance_constraints

"""
$(TYPEDSIGNATURES)

Run the Community Flux Balance Analysis on several models. All arguments are
forwarded to [`community_flux_balance_constraints`](@ref) which constructs the
model; this function returns the solved model.
"""
community_flux_balance_analysis(args...; kwargs...) = frontend_optimized_values(
    community_flux_balance_constraints,
    args...;
    objective = x -> x.community_biomass.value,
    kwargs...,
)

export community_flux_balance_analysis

# TODO this might eventually deserve the QP community-FBA version

"""
$(TYPEDSIGNATURES)

Build a community of `models` of fixed total biomass `growth`, where the
community abundances are represented as variables. This allows finding feasible
community compositions at a given growth rate, following the methodology of
SteadyCom algorithm.

`models` is an interable of pairs in the shape of `(name, model)`, where a
`model` can be a constraint tree or an abstract metabolic model. The models are
connected by interfaces and bounded just like in
`community_flux_balance_constraints` (as the main difference, this function
does not accept abundances, but requires knowledge of growth in advance).

The output constraint tree contains the target growth value as `growth`,
relative abundances in subtree `abundances`, all values of the original
models in subtree `community`, and SteadyCom-style diluted constraints in
subtree `diluted_constraints`.
"""
function community_composition_balance_constraints(
    models,
    growth,
    community_exchange_bounds = Dict();
    interface = :identifier_prefixes,
    interface_exchanges = x -> x.interface.exchanges,
    interface_biomass = x -> x.interface.biomass,
    default_community_exchange_bound = nothing,
)
    @assert length(models) >= 1 "at least one community member is required"
    @assert growth >= 0 "growth must not be negative"

    bounds_lookup = Dict(community_exchange_bounds)

    make_member(k, m) = throw(DomainError(m, "unsupported community member type ($k)"))
    make_member(k, m::A.AbstractFBCModel) =
        make_member(k, flux_balance_constraints(m; interface))
    make_member(_, c::C.ConstraintTree) = (c, interface_exchanges(c))

    constraints = interface_constraints(
        (Symbol(k) => make_member(k, m) for (k, m) in models);
        out_interface = :community_exchanges,
        out_balance = :community_balance,
        wrap = x -> :community^x,
        bound = x ->
            get(bounds_lookup, String(last(x)), default_community_exchange_bound),
    )

    constraints +=
        :abundances^C.variables(
            keys = collect(keys(constraints.community)),
            bounds = C.Between(0, 1),
        )

    return C.ConstraintTree(
        (kv for kv in constraints if first(kv) != :community)...,
        :community => remove_bounds(constraints.community),
        :diluted_constraints => C.ConstraintTree(
            k => value_scaled_bound_constraints(v, constraints.abundances[k].value) for
            (k, v) in constraints.community
        ),
        :total_abundance_constraint => C.Constraint(sum_value(constraints.abundances), 1.0),
        :growth => C.Constraint(C.LinearValue(growth)),
        :biomass_constraints => C.ConstraintTree(
            k => equal_value_constraint(
                sum_value(interface_biomass(constraints.community[k])),
                constraints.abundances[k].value * growth,
            ) for (k, v) in constraints.community
        ),
    )
end

export community_composition_balance_constraints

"""
$(TYPEDSIGNATURES)

Build a variable-composition community from given models, and scan for the
maximum achievable growth. The scanning proceeds by halving the
interval `(0, maximum_growth)`, until the interval size becomes smaller than
the desired `tolerance`.

Positional and keyword arguments are forwarded to
[`community_composition_balance_constraints`](@ref). `optimizer`, `settings`,
`tolerance` and `output` are forwarded to [`feasibility_threshold`](@ref).

Function `output` may be specified to restrict the reported result: For
example, to just return the best achieved growth and abundances, one may use
`output = x -> :growth^x.growth * :abundances^x.abundances`.
"""
function community_composition_balance_analysis(
    models,
    maximum_growth,
    args...;
    tolerance = 1e-4,
    optimizer,
    settings = [],
    output = identity,
    kwargs...,
)
    @assert maximum_growth >= 0 "maximum_growth must not be negative"

    return feasibility_threshold(
        0.0,
        Float64(maximum_growth);
        tolerance,
        optimizer,
        settings,
        output,
    ) do growth
        community_composition_balance_constraints(models, growth, args...; kwargs...)
    end
end

export community_composition_balance_analysis

"""
$(TYPEDSIGNATURES)

Run a SteadyCom-style analysis, and return the variability of the community
composition at the given growth rate.

Positional and keyword arguments are forwarded to
[`community_composition_balance_constraints`](@ref); `optimizer`, `settings`
and `workers` are internally forwarded to [`constraints_variability`](@ref).
"""
function community_composition_variability_analysis(
    models,
    growth,
    args...;
    optimizer,
    settings = [],
    workers = D.workers(),
    kwargs...,
)
    constraints =
        community_composition_balance_constraints(models, growth, args...; kwargs...)
    return constraints_variability(
        constraints,
        constraints.abundances;
        optimizer,
        workers,
        settings,
    )
end

export community_composition_variability_analysis
