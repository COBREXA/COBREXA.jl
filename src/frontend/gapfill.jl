
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

"""
$(TYPEDSIGNATURES)

Make a gap-filling system from a FBC model with gaps and an universal
FBC model that contains reactions to be added into the original model.

The output system will be constrainted to reach at least `objective_target`
flux through the objective function. Generally, this should be set to an
arbitrary small value such as `0.05`.

`universal_reaction_cost` should assign a numeric cost of inclusion of each of
the reactions in the `universal_model`; by default all are assigned equal
weight of `1`. `max_cost` puts an optional maximum limit on the cost, which may
help the solver to avoid exploring unnecessarily complex solutions.
`known_fills` may contain previous solutions of the same system; these will be
made infeasible in the output constraint system in order to allow discovery of
other ones.

Additional arguments are forwarded to `flux_balance_constraints` that converts
`model` to constraints.
"""
function gap_filling_constraints(
    model::A.AbstractFBCModel,
    universal_model::A.AbstractFBCModel,
    objective_target::Float64;
    universal_reaction_cost = _ -> 1.0,
    max_cost = Inf,
    known_fills::Vector{C.Tree{Float64}} = C.Tree{Float64}[],
    kwargs...,
)
    m = flux_balance_constraints(model; kwargs...)
    m.objective.bound = C.Between(objective_target, Inf)
    s = m.flux_stoichiometry
    delete!(m, :flux_stoichiometry)
    u = flux_balance_constraints(universal_model)
    gap_filling_constraints(;
        system = m,
        stoichiometry = s,
        universal_fluxes = u.fluxes,
        universal_stoichiometry = u.flux_stoichiometry,
        flux_cost = i -> universal_reaction_cost(String(i)),
        max_cost,
        known_fills,
    )
end

"""
$(TYPEDSIGNATURES)

Make a gap-fillign system from a non-solving `system` with a separated
`stoichiometry`, filling in possible fluxes from `universal_fluxes` that are
balanced with `universal_stoichiometry`

`flux_cost` can be used to assign a weight to each given universal flux; the
total cost is bounded by `max_cost`.

`known_fills` may contain previous solutions to be avoided; these are processed
by [`gap_filling_known_fill_constraint`](@ref) and attached to the system.

`stoichiometry` needs to be extended to construct the final constraints, thus
it should not be present in the original `system`.
"""
function gap_filling_constraints(;
    system::C.ConstraintTree,
    stoichiometry::C.ConstraintTree,
    universal_fluxes::C.ConstraintTree,
    universal_stoichiometry::C.ConstraintTree,
    flux_cost = _ -> 1.0,
    max_cost = Inf,
    known_fills::Vector{C.Tree{Float64}} = C.Tree{Float64}[],
)
    joined =
        C.ConstraintTree(:system => system, :stoichiometry => stoichiometry) +
        :universal^C.ConstraintTree(
            :fluxes => universal_fluxes,
            :stoichiometry => universal_stoichiometry,
        ) +
        :fill_flags^C.variables_for(universal_fluxes) do _
            Switch(0, 1)
        end

    return C.ConstraintTree(
        :system => joined.system,
        :universal_fluxes => joined.universal.fluxes,
        :universal_flux_bounds => C.zip(joined.universal.fluxes, joined.fill_flags) do x, b
            if x.bound isa C.Between
                C.ConstraintTree(
                    :lower => C.Constraint(
                        x.value - x.bound.lower * b.value,
                        C.Between(0, Inf),
                    ),
                    :upper => C.Constraint(
                        x.value - x.bound.upper * b.value,
                        C.Between(-Inf, 0),
                    ),
                )
            elseif x.bound isa C.EqualTo
                C.Constraint(x.value - x.bound.equal_to * b.value, 0)
            elseif isnothing(x.bound)
                C.ConstraintTree()
            else
                throw(DomainError(x.bound, "unsupported flux bound"))
            end
        end,
        :stoichiometry =>
            C.merge(joined.stoichiometry, joined.universal.stoichiometry) do a, b
                ismissing(a) && return b
                ismissing(b) && return a
                @assert a.bound isa C.EqualTo &&
                        a.bound.equal_to == 0.0 &&
                        b.bound isa C.EqualTo &&
                        b.bound.equal_to == 0.0 "Stoichiometries in both systems must only contain equal-to-zero-bounded constraints"
                C.Constraint(a.value + b.value, 0.0)
            end,
        :fill_flags => joined.fill_flags,
        (
            Symbol(:known_fills_, i) =>
                gap_filling_known_fill_constraint(joined.fill_flags, kf) for
            (i, kf) in enumerate(known_fills)
        )...,
        :n_filled => C.Constraint(
            C.sum(
                (flux_cost(k) * v.value for (k, v) in joined.fill_flags),
                init = zero(C.LinearValue),
            ),
            C.Between(0, max_cost),
        ),
    )
end

export gap_filling_constraints

"""
$(TYPEDSIGNATURES)

Produce a constraint that can be added to the system made by
[`gap_filling_constraints`](@ref) to avoid repeating of a solution that
includes reactions already generated by another solution, as represented by the
solved `fill_flags`.

Parameter `fill_flags` are the gapfilling flags of the given constraint system,
parameter `known_flags` is expected to contain the solved `fill_flags` for the
solution that is to be avoided.
"""
gap_filling_known_fill_constraint(
    fill_flags::C.ConstraintTree,
    known_flags::C.Tree{Float64},
) = C.Constraint(
    C.sum(
        values(C.zip(fill_flags, known_flags, C.Value) do f, k
            k - f.value * k
        end),
        init = zero(C.LinearValue),
    ),
    (1 - eps(), Inf),
)

export gap_filling_known_fill_constraint

"""
$(TYPEDSIGNATURES)

Run the gap-filling analysis on a constraint system specified by
[`gap_filling_constraints`](@ref).
"""
gap_filling_analysis(args...; kwargs...) = frontend_optimized_values(
    gap_filling_constraints,
    args...;
    objective = x -> x.n_filled.value,
    sense = Minimal,
    kwargs...,
)

export gap_filling_analysis
