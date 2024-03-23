
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

# TODO this might eventually deserve the QP version

# TODO steadycom

"""
$(TYPEDSIGNATURES)

TODO
"""
function community_flux_balance_constraints(
    model_abundances,
    community_exchange_bounds;
    interface = :identifier_prefixes,
    interface_exchanges = x -> x.interface.exchanges,
    interface_biomass = x -> x.interface.biomass,
)
    @assert length(model_abundances) >= 1 "at least one community member is required"
    # TODO test if abundances sum to 1?

    bounds_lookup = Dict(community_exchange_bounds)

    constraints = interface_constraints(
        (
            Symbol(k) => let
                c = flux_balance_constraints(m; interface)
                (c, interface_exchanges(c), a)
            end for (k, (m, a)) in model_abundances
        );
        out_interface = :community_exchanges,
        out_balance = :community_balance,
        bound = x -> get(bounds_lookup, String(last(x)), nothing),
    )

    growth_sums = [
        Symbol(k) => C.Constraint(sum_value(interface_biomass(constraints[Symbol(k)])))
        for (k, _) in model_abundances
    ]

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

TODO
"""
community_flux_balance_analysis(args...; kwargs...) = frontend_optimized_values(
    community_flux_balance_constraints,
    args...;
    objective = x -> x.community_biomass.value,
    kwargs...,
)

export community_flux_balance_analysis
