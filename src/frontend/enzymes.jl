
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
$(TYPEDEF)

A simple struct storing information about the isozyme composition, including
subunit stoichiometry and turnover numbers. Use with
[`enzyme_constrained_flux_balance_analysis`](@ref).

If either of the turnover numbers is `nothing`, no bounds are added for that
direction; i.e., the reaction is assumed not to require enzymes to proceed in
that direction. If you want to disable the reaction in the given direction
instead, either use `0` as a turnover number, or better put bounds directly on
reaction flux.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct IsozymeT{T<:Real}
    """
    Mapping of gene product identifiers ("genes" in FBC model nomenclature)
    to their relative amount required to construct one unit of the isozyme.
    """
    gene_product_stoichiometry::Dict{String,T}

    "Turnover number for this isozyme catalyzing the forward direction of the
    reaction, or `nothing` if isozyme is not needed."
    kcat_forward::Maybe{T} = nothing

    "Turnover number for this isozyme catalyzing the reverse direction of the
    reaction, or `nothing` if isozyme is not needed."
    kcat_reverse::Maybe{T} = nothing
end

export IsozymeT

"""
$(TYPEDEF)

Shortcut for `[IsozymeT](@ref){Float64}`.
"""
const Isozyme = IsozymeT{Float64}

export Isozyme


"""
$(TYPEDSIGNATURES)

Expand the `capacity` argument as given to
[`enzyme_constrained_flux_balance_constraints`](@ref) and
[`simplified_enzyme_constrained_flux_balance_constraints`](@ref) into a form
accepted by [`enzyme_constraints`](@ref) and
[`simplified_enzyme_constraints`](@ref) (respectively).

By default, `Bound`s are kept intact, `Real` values are converted to a fixed
interval between a zero and the value. All other values are assumed to be lists
of capacities. (See [`expand_enzyme_capacity_bound`](@ref) for translation of
actual bounds).

Overloading this function (or [`expand_enzyme_capacity_bound`](@ref)) gives a
way to simplify the interface of the functions by accomodating custom capacity
types.

The second argument is provided to this function as a list of full scope of the
capacities it can work with, by default "all capacities".
"""
expand_enzyme_capacity(x, all) =
    return [:total_capacity => (all, expand_enzyme_capacity_bound(x))]

"""
$(TYPEDSIGNATURES)

Overload of [`expand_enzyme_capacity`](@ref) for all `Dict`-like iterables.
"""
expand_enzyme_capacity(x::Union{Vector{Pair},Dict}, _) =
    return expand_enzyme_capacity_iterable(x)

"""
$(TYPEDSIGNATURES)

Overload of [`expand_enzyme_capacity`](@ref) that provides compatibility with
the earlier capacity specifications (using triples instead of pairs).
"""
expand_enzyme_capacity(x::Vector{<:Tuple}, _) =
    return expand_enzyme_capacity_iterable(id => (grp, cap) for (id, grp, cap) in x)

export expand_enzyme_capacity

"""
$(TYPEDSIGNATURES)

Internal helper for implementation of [`expand_enzyme_capacity`](@ref) over
iterable `Dict`-like objects.
"""
expand_enzyme_capacity_iterable(x) = return [
    Symbol(id) => (grp, expand_enzyme_capacity_bound(cap)) for (id, (grp, cap)) in x
]

"""
$(TYPEDSIGNATURES)

Expand a single capacity bound for use in enzyme-constrained models.
Overloading this function provides additional ways to interpret the capacity
specifications. Typically used via [`expand_enzyme_capacity`](@ref).
"""
expand_enzyme_capacity_bound(x::Real) = return (zero(x), x)

"""
$(TYPEDSIGNATURES)

By default, [`expand_enzyme_capacity_bound`](@ref) leaves all `Bound`s intact.
This overload implements this property.
"""
expand_enzyme_capacity_bound(x::C.Bound) = return x

export expand_enzyme_capacity_bound

"""
$(TYPEDSIGNATURES)

Construct a enzyme-constrained flux-balance constraint system,
following the method in GECKO algorithm (refer to: *Sánchez, Benjamín J., et
al. "Improving the phenotype predictions of a yeast genome‐scale metabolic
model by incorporating enzymatic constraints." Molecular systems biology 13.8
(2017): 935*).

The enzyme mass constraints depend primarily on the available *isozymes*, given
in parameter `reaction_isozymes`, which is a mapping of reaction identifiers to
descriptions of [`Isozyme`](@ref)s that may catalyze the particular reactions.
The isozymes are built from gene products, the mass of which is specified by
`gene_product_molar_masses`. In total, the amount of gene product building
material is limited by `capacity`.

`capacity` may be a single number, which sets the mass limit for "all described
enzymes". Alternatively, `capacity` may be a vector of identifier-genes-limit
triples that together form a constraint (identified by the given identifier)
that limits the total sum of the listed genes to the given limit. The
interpretation of `capacity` is implemented (and can be extended) via
[`expand_enzyme_capacity`](@ref).

`interface` and `interface_name` are forwarded to
[`flux_balance_constraints`](@ref).
"""
function enzyme_constrained_flux_balance_constraints(
    model::A.AbstractFBCModel;
    reaction_isozymes::Dict{String,Dict{String,IsozymeT{R}}},
    gene_product_molar_masses::Dict{String,Float64},
    capacity,
    interface::Maybe{Symbol} = nothing,
    interface_name = :interface,
) where {R<:Real}
    # prepare some accessor functions for the later stuff
    # TODO: might be nicer to somehow parametrize the fwd/rev directions out.
    # Also there is a lot of conversion between symbols and strings, might be
    # nicer to have that sorted out in some better way.
    function isozyme_forward_ids(rid)
        haskey(reaction_isozymes, String(rid)) || return nothing
        return [
            Symbol(k) for
            (k, i) in reaction_isozymes[String(rid)] if !isnothing(i.kcat_forward)
        ]
    end
    function isozyme_reverse_ids(rid)
        haskey(reaction_isozymes, String(rid)) || return nothing
        return [
            Symbol(k) for
            (k, i) in reaction_isozymes[String(rid)] if !isnothing(i.kcat_reverse)
        ]
    end
    kcat_forward(rid, iso_id) = reaction_isozymes[String(rid)][String(iso_id)].kcat_forward
    kcat_reverse(rid, iso_id) = reaction_isozymes[String(rid)][String(iso_id)].kcat_reverse
    isozyme_gene_product_stoichiometry(rid, iso_id) = Dict(
        Symbol(k) => v for (k, v) in
        reaction_isozymes[String(rid)][String(iso_id)].gene_product_stoichiometry
    )
    gene_ids = Symbol.(keys(gene_product_molar_masses))
    gene_product_molar_mass(gid) = get(gene_product_molar_masses, String(gid), 0.0)

    # allocate all variables and build the system
    constraints = flux_balance_constraints(model; interface, interface_name)

    constraints +=
        :fluxes_forward^unsigned_positive_contribution_variables(constraints.fluxes)
    constraints *=
        :fluxes_reverse^unsigned_negative_contribution_constraints(
            constraints.fluxes,
            constraints.fluxes_forward,
        )

    constraints += enzyme_variables(;
        fluxes_forward = constraints.fluxes_forward,
        fluxes_reverse = constraints.fluxes_reverse,
        isozyme_forward_ids,
        isozyme_reverse_ids,
    )

    return constraints * enzyme_constraints(;
        fluxes_forward = constraints.fluxes_forward,
        fluxes_reverse = constraints.fluxes_reverse,
        isozyme_forward_amounts = constraints.isozyme_forward_amounts,
        isozyme_reverse_amounts = constraints.isozyme_reverse_amounts,
        kcat_forward,
        kcat_reverse,
        isozyme_gene_product_stoichiometry,
        gene_product_molar_mass,
        capacity_limits = expand_enzyme_capacity(capacity, gene_ids),
    )
end

export enzyme_constrained_flux_balance_constraints

"""
$(TYPEDSIGNATURES)

Perform the enzyme-constrained flux balance analysis on the `model` and return
the solved constraint system.

Arguments are forwarded to
[`enzyme_constrained_flux_balance_constraints`](@ref); solver configuration
arguments are forwarded to [`optimized_values`](@ref).
"""
enzyme_constrained_flux_balance_analysis(model::A.AbstractFBCModel; kwargs...) =
    frontend_optimized_values(
        enzyme_constrained_flux_balance_constraints,
        model;
        objective = x -> x.objective.value,
        kwargs...,
    )

export enzyme_constrained_flux_balance_analysis

"""
$(TYPEDSIGNATURES)

Like [`enzyme_constrained_flux_balance_constraints`](@ref), but automatically
selects a single "fastest" isozyme for each reaction direction. In turn, the
system requires much less variables in the constraint system description, and
usually solves more efficiently, for the price of possibly finding suboptimal
solutions. The method follows the SMOMENT algorithm described in *Bekiaris,
P.S., Klamt, S. Automatic construction of metabolic models with enzyme
constraints. BMC Bioinformatics 21, 19 (2020).
https://doi.org/10.1186/s12859-019-3329-9*.

Arguments are as with [`enzyme_constrained_flux_balance_constraints`](@ref),
with a major difference in `capacity` handling: the identifier lists contain
reactions identifiers (i.e., keys of `fluxes` in the constraint tree), instead
of the gene product identifiers.
"""
function simplified_enzyme_constrained_flux_balance_constraints(
    model;
    reaction_isozymes::Dict{String,Dict{String,IsozymeT{R}}},
    gene_product_molar_masses::Dict{String,Float64},
    capacity,
    interface::Maybe{Symbol} = nothing,
    interface_name = :interface,
) where {R<:Real}
    # TODO this deserves a rewrite once more -- Isozyme struct is hiding a bit
    # too much of uncertainty for the code of this thing to be short and
    # concise...maybe we should have an isozyme with only one kcat which is
    # never missing so that the nothing-handling code and duplicate getters can
    # disappear?

    # prepare getters and data

    isozyme_mass(i) = C.sum(
        (
            (
                haskey(gene_product_molar_masses, gid) ?
                coeff * gene_product_molar_masses[gid] :
                throw(DomainError(gid, "missing a required gene product molar mass"))
            ) for (gid, coeff) in i.gene_product_stoichiometry
        ),
        init = 0.0,
    )

    min_isozyme_cost_forward = Dict(
        Symbol(rid) => argmin(
            last,
            (i, isozyme_mass(i) / i.kcat_forward) for
            (_, i) in is if !isnothing(i.kcat_forward)
        ) for (rid, is) in reaction_isozymes if
        any(!isnothing(i.kcat_forward) for (_, i) in is)
    )
    min_isozyme_cost_reverse = Dict(
        Symbol(rid) => argmin(
            last,
            (i, isozyme_mass(i) / i.kcat_reverse) for
            (_, i) in is if !isnothing(i.kcat_reverse)
        ) for (rid, is) in reaction_isozymes if
        any(!isnothing(i.kcat_forward) for (_, i) in is)
    )

    # allocate the model and variables

    constraints = flux_balance_constraints(model; interface, interface_name)

    constraints +=
        :fluxes_forward^unsigned_positive_contribution_variables(constraints.fluxes)
    constraints *=
        :fluxes_reverse^unsigned_negative_contribution_constraints(
            constraints.fluxes,
            constraints.fluxes_forward,
        )

    # connect everything with constraints

    return constraints *
           :capacity_limits^simplified_enzyme_constraints(;
               fluxes_forward = constraints.fluxes_forward,
               fluxes_reverse = constraints.fluxes_reverse,
               mass_cost_forward = rid ->
                   maybemap(last, get(min_isozyme_cost_forward, rid, nothing)),
               mass_cost_reverse = rid ->
                   maybemap(last, get(min_isozyme_cost_reverse, rid, nothing)),
               capacity_limits = expand_enzyme_capacity(capacity, keys(constraints.fluxes)),
           ) *
           :gene_product_amounts^simplified_isozyme_gene_product_amount_constraints(
               (
                   constraints.fluxes_forward,
                   rid -> maybemap(
                       ic -> (
                           Dict(
                               Symbol(k) => v for
                               (k, v) in first(ic).gene_product_stoichiometry
                           ),
                           first(ic).kcat_forward,
                       ),
                       get(min_isozyme_cost_forward, rid, nothing),
                   ),
               ),
               (
                   constraints.fluxes_reverse,
                   rid -> maybemap(
                       ic -> (
                           Dict(
                               Symbol(k) => v for
                               (k, v) in first(ic).gene_product_stoichiometry
                           ),
                           first(ic).kcat_reverse,
                       ),
                       get(min_isozyme_cost_reverse, rid, nothing),
                   ),
               ),
           )
    # TODO add some generic functionality that allows people to add additional
    # capacity group bounds over gene product amounts. That code can be shared
    # with enzyme_constraints.
end

export simplified_enzyme_constrained_flux_balance_constraints

"""
$(TYPEDSIGNATURES)

Perform the enzyme-constrained flux balance analysis on the `model` and return
the solved constraint system.

Arguments are forwarded to
[`simplified_enzyme_constrained_flux_balance_constraints`](@ref); solver
configuration arguments are forwarded to [`optimized_values`](@ref).
"""
simplified_enzyme_constrained_flux_balance_analysis(model::A.AbstractFBCModel; kwargs...) =
    frontend_optimized_values(
        simplified_enzyme_constrained_flux_balance_constraints,
        model;
        objective = x -> x.objective.value,
        kwargs...,
    )

export simplified_enzyme_constrained_flux_balance_analysis
