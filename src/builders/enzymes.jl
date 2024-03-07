
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

Create a `ConstraintTree` with variables for isozyme contributions to reaction
fluxes. The tree has 2 levels: the first contains all reaction flux IDs that
have isozymes, the second contains the isozyme IDs for each reaction flux.

`fluxes` should be anything that can be iterated to give reaction flux IDs.

`flux_isozymes` is a function that, for a given reaction flux ID, returns
anything iterable that contains the isozyme IDs for the given reaction flux.
Returning an empty iterable prevents allocating the subtree for the given flux.
"""
isozyme_amount_variables(fluxes, flux_isozymes) = sum(
    (
        f^C.variables(keys = fis, bounds = C.Between(0, Inf)) for
        (f, fis) in ((f, flux_isozymes(f)) for f in fluxes) if
        !isnothing(fis) && !isempty(fis)
    ),
    init = C.ConstraintTree(),
)

export isozyme_amount_variables

"""
$(TYPEDSIGNATURES)

A constraint tree that sums up partial contributions of reaction isozymes to
the fluxes of reactions.

For practical purposes, both fluxes and isozymes are here considered to be
unidirectional, i.e., one would typically apply this twice to constraint both
"forward" and "reverse" fluxes.

Function `kcat` should return the kcat value for a given reaction and isozyme
(IDs of which respectively form the 2 parameters for each call).
"""
function isozyme_flux_constraints(
    isozyme_amounts::C.ConstraintTree,
    fluxes::C.ConstraintTree,
    kcat,
)
    C.ConstraintTree(
        rid => C.Constraint(
            sum(kcat(rid, iid) * i.value for (iid, i) in ri if !isnothing(kcat(rid, iid))) - fluxes[rid].value,
            0.0,
        ) for (rid, ri) in isozyme_amounts if haskey(fluxes, rid)
    )
end

export isozyme_flux_constraints

"""
$(TYPEDSIGNATURES)

A constraint tree that computes the gene product amounts from given isozyme
amounts their multiplicities (aka. stoichiometries, protein units, ...) given
by `isozyme_stoichiometry`.

Values in ConstraintTree `gene_product_amounts` should describe the gene
product allocations.  Allocation for the isozyme is ignored if the gene product
is missing in `gene_product_amounts`.

`isozyme_amounts` is an iterable that contains several `ConstraintTree`s that
describe the allocated isozyme amounts (typically these would be created by
[`isozyme_amount_variables`](@ref). The multiple trees may describe several
different kinds of isozyme use, e.g., you can use it to pass in both forward-
and reverse-direction amounts at once. To only use a single tree, use an
uni-tuple: `isozyme_amounts = tuple(my_tree)`.

Parameter function `isozyme_stoichiometry` gets called with a reaction and
isozyme IDs as given by the isozyme amount trees. It should return `nothing` in
case there's no information -- in such case, the isozyme is *not going to be
included* in the calculation of gene product mass.
"""
function isozyme_gene_product_amount_constraints(isozyme_amounts, isozyme_stoichiometry)
    res = C.ConstraintTree()
    # This needs to "invert" the stoichiometry mapping,
    # so we patch up a fresh new CT in place.
    for iss in isozyme_amounts
        for (rid, is) in iss
            for (iid, i) in is
                gpstoi = isozyme_stoichiometry(rid, iid)
                isnothing(gpstoi) && continue
                for (gp, stoi) in gpstoi
                    if haskey(res, gp)
                        res[gp].value += i.value * stoi
                    else
                        res[gp] = C.Constraint(i.value * stoi)
                    end
                end
            end
        end
    end
    res
end

export isozyme_gene_product_amount_constraints

"""
$(TYPEDSIGNATURES)

Returns a constraint tree with enzyme capacity constraints, added for reactions
in `fluxes_forward` and `fluxes_reverse`. This is used to construct the
constraint system in [`enzyme_constrained_flux_balance_constraints`](@ref).

Parameter function `isozyme_ids` takes a reaction ID and returns `nothing` if
the reaction does not have isozymes associated with it, or an iterable
container of all the isozyme IDs for that reaction (as `Symbol`s).

Parameters `isozyme_forward_ids` and `isozyme_reverse_ids` can be used to
fine-tune the generated isozymes in either direction; both default to
`isozyme_ids`.

The keys in the output constraint tree can be customized by setting
`isozyme_forward_amounts_name`, `isozyme_reverse_amounts_name` and
`gene_product_amounts_name`.
"""
enzyme_variables(;
    fluxes_forward::C.ConstraintTree,
    fluxes_reverse::C.ConstraintTree,
    isozyme_ids = _ -> nothing,
    isozyme_forward_ids = isozyme_ids,
    isozyme_reverse_ids = isozyme_ids,
    isozyme_forward_amounts_name = :isozyme_forward_amounts,
    isozyme_reverse_amounts_name = :isozyme_reverse_amounts,
) =
    isozyme_forward_amounts_name^isozyme_amount_variables(
        Symbol.(keys(fluxes_forward)),
        isozyme_forward_ids,
    ) +
    isozyme_reverse_amounts_name^isozyme_amount_variables(
        Symbol.(keys(fluxes_reverse)),
        isozyme_reverse_ids,
    )

export enzyme_variables

"""
$(TYPEDSIGNATURES)

Connect variables returned by [`enzyme_variables`](@ref) to unidirectional
fluxes. This is used to construct the contraint system for
[`enzyme_constrained_flux_balance_constraints`](@ref).

Parameters `fluxes_forward`, `fluxes_reverse`, `isozyme_forward_amounts`,
`isozyme_reverse_amounts` and `gene_product_amounts` should correspond to
parameters and results of [`enzyme_variables`](@ref).

Further, parameter functions `kcat_forward` and `kcat_reverse` specify the
turnover numbers for reaction and isozyme IDs given in parameters;
`isozyme_gene_product_stoichiometry` specifies the composition of the
reaction-isozyme IDs given in parameter by returning an interable mapping of
gene product IDs to numbers (such as `Dict{Symbol, Float64}`), and
`gene_product_molar_mass` specifies a numeric mass for a given gene product ID.
All parameter functions may return `nothing`, at which point the given
object is considered nonexistent and is omitted from constraints.

`capacity_limits` is an interable container of triples `(limit_id,
gene_product_ids, capacity_bound)` which are converted to a constraint
identified by the `limit_id` that limits the total mass of `gene_product_ids`
(which is any iterable container) by `capacity_bound`.
"""
function enzyme_constraints(;
    fluxes_forward::C.ConstraintTree,
    fluxes_reverse::C.ConstraintTree,
    isozyme_forward_amounts::C.ConstraintTree,
    isozyme_reverse_amounts::C.ConstraintTree,
    kcat_forward = (_, _) -> nothing,
    kcat_reverse = (_, _) -> nothing,
    isozyme_gene_product_stoichiometry = (_, _) -> nothing,
    gene_product_molar_mass = _ -> nothing,
    capacity_limits = [],
    isozyme_flux_forward_balance_name = :isozyme_flux_forward_balance,
    isozyme_flux_reverse_balance_name = :isozyme_flux_reverse_balance,
    gene_product_amounts_name = :gene_product_amounts,
    gene_product_capacity_name = :gene_product_capacity,
)
    gene_product_amounts = isozyme_gene_product_amount_constraints(
        (isozyme_forward_amounts, isozyme_reverse_amounts),
        isozyme_gene_product_stoichiometry,
    )

    return isozyme_flux_forward_balance_name^isozyme_flux_constraints(
               isozyme_forward_amounts,
               fluxes_forward,
               kcat_forward,
           ) *
           isozyme_flux_reverse_balance_name^isozyme_flux_constraints(
               isozyme_reverse_amounts,
               fluxes_reverse,
               kcat_reverse,
           ) *
           gene_product_amounts_name^gene_product_amounts *
           gene_product_capacity_name^C.ConstraintTree(
               id => C.Constraint(;
                   value = sum(
                       gene_product_amounts[gp].value * gpmm for
                       (gp, gpmm) in ((gp, gene_product_molar_mass(gp)) for gp in gps) if
                       !isnothing(gpmm) && haskey(gene_product_amounts, gp)
                   ),
                   bound,
               ) for (id, gps, bound) in capacity_limits
           )
end

export enzyme_constraints

"""
$(TYPEDSIGNATURES)

Like [`isozyme_gene_product_amount_constraints`](@ref), but works with the
"simplified" view where each reaction has an uniquely determined catalytic
isozyme, as with [`simplified_enzyme_constraints`](@ref).

As the main difference, the arguments are tuples that contain first the
constraint tree without the "isozyme" layer (i.e., fluxes), and second a
function that returns the gene product stoichiometry and the turnover number
(again in a tuple) for the given flux identifier.
"""
function simplified_isozyme_gene_product_amount_constraints(x...)
    res = C.ConstraintTree()
    for (iss, stoikcatfn) in x
        for (rid, i) in iss
            x = stoikcatfn(rid)
            isnothing(x) && continue
            (gpstoi, kcat) = x
            isnothing(kcat) && continue
            for (gp, stoi) in gpstoi
                if haskey(res, gp)
                    res[gp].value += i.value * stoi / kcat
                else
                    res[gp] = C.Constraint(i.value * stoi / kcat)
                end
            end
        end
    end
    res
end

export simplified_isozyme_gene_product_amount_constraints

"""
$(TYPEDSIGNATURES)

Build a constraint system that bounds fluxes according to their enzyme mass
requirements, with respect to per-reaction enzyme mass costs.

Parameter functions `mass_cost_forward` and `mass_cost_reverse` take a flux ID
(corresponding to a flux in `fluxes_forward` and `fluxes_reverse`) and return
the enzyme mass required to catalyze one "unit" of reaction in the forward or
reverse direction, respectively. Returning `nothing` ignores the mass cost.

`capacity_limits` is an iterable container of triples `(limit_id, flux_ids,
bound)`, which creates the capacity bounds over groups of fluxes (in the same
manner as for gene products in [`enzyme_constraints`](@ref)).
"""
function simplified_enzyme_constraints(;
    fluxes_forward::C.ConstraintTree,
    fluxes_reverse::C.ConstraintTree,
    mass_cost_forward = _ -> nothing,
    mass_cost_reverse = _ -> nothing,
    capacity_limits = [],
)
    function contribution(fl, cost, id)
        c = cost(id)
        (isnothing(c) || !haskey(fl, id)) && return zero(C.LinearValue)
        return c * fl[id].value
    end

    return C.ConstraintTree(
        id => C.Constraint(;
            value = sum(
                contribution(fluxes_forward, mass_cost_forward, f) for f in fs;
                init = zero(C.LinearValue),
            ) + sum(
                contribution(fluxes_reverse, mass_cost_reverse, f) for f in fs;
                init = zero(C.LinearValue),
            ),
            bound,
        ) for (id, fs, bound) in capacity_limits
    )
end

export simplified_enzyme_constraints
