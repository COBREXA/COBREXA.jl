
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

A constraint tree that models the content of the given instance of
`AbstractFBCModel`.

The constructed tree contains subtrees `fluxes` (with the reaction-defining
"variables") and `flux_stoichiometry` (with the metabolite-balance-defining
constraints), and a single constraint `objective` thad describes the objective
function of the model.

Optionally if `interface` is specified, an "interface block" will be created
within the constraint tree for later use as a "module" in creating bigger
models (such as communities) using [`interface_constraints`](@ref). The
possible parameter values include:
- `nothing` -- default, no interface is created
- `:sbo` -- the interface gets created from model's SBO annotations)
- `:identifier_prefixes` -- the interface is guesstimated from commonly
  occurring adhoc reaction ID prefixes used in contemporary models
- `:boundary` -- the interface is created from all reactions that either only
  consume or only produce metabolites

Output interface name can be set via `interface_name`.

See [`Configuration`](@ref) for fine-tuning the default interface creation.
"""
function flux_balance_constraints(
    model::A.AbstractFBCModel;
    interface::Maybe{Symbol} = nothing,
    interface_name = :interface,
)
    rxn_strings = A.reactions(model)
    rxns = Symbol.(rxn_strings)
    mets = Symbol.(A.metabolites(model))
    lbs, ubs = A.bounds(model)
    stoi = A.stoichiometry(model)
    bal = A.balance(model)
    obj = A.objective(model)
    cpls = Symbol.(A.couplings(model))
    cT = SparseArrays.SparseMatrixCSC(A.coupling(model)')
    clbs, cubs = A.coupling_bounds(model)

    # The iteration through stoichiometry would be better done with eachrow(),
    # unfortunately it seems to be enormously inefficient on column-major
    # matrix formats, and findnz() on eachrow&eachcol output isn't very
    # backwards compatible (at least not to 1.6). We thus enforce having a
    # SparseMatrixCSC here and walk it manually.
    stoiT = SparseArrays.SparseMatrixCSC(stoi')

    constraints = C.ConstraintTree(
        :fluxes^C.variables(keys = rxns, bounds = zip(lbs, ubs)) *
        :flux_stoichiometry^C.ConstraintTree(
            met => C.Constraint(
                value = let i = stoiT.colptr[row_idx], e = stoiT.colptr[row_idx+1] - 1

                    C.LinearValue(idxs = stoiT.rowval[i:e], weights = stoiT.nzval[i:e])
                end,
                bound = C.EqualTo(b),
            ) for (met, row_idx, b) in zip(mets, 1:stoiT.n, bal)
        ) *
        :coupling^C.ConstraintTree(
            cpl => C.Constraint(
                value = let i = cT.colptr[row_idx], e = cT.colptr[row_idx+1] - 1
                    C.LinearValue(idxs = cT.rowval[i:e], weights = cT.nzval[i:e])
                end,
                bound = clb == cub ? C.EqualTo(clb) : C.Between(clb, cub),
            ) for (cpl, row_idx, clb, cub) in zip(cpls, 1:cT.n, clbs, cubs)
        ) *
        :objective^C.Constraint(C.LinearValue(SparseArrays.sparse(obj))),
    )

    add_interface(sym, flt) =
        any(flt) && (
            constraints *=
                interface_name^sym^C.ConstraintTree(
                    r => constraints.fluxes[r] for r in rxns[flt]
                )
        )
    if interface == :sbo
        sbod(sbos, rid) = any(in(sbos), get(A.reaction_annotations(model, rid), "sbo", []))
        add_interface(:exchanges, sbod.(Ref(configuration.exchange_sbos), rxn_strings))
        add_interface(:biomass, sbod.(Ref(configuration.biomass_sbos), rxn_strings))
        add_interface(
            :atp_maintenance,
            sbod.(Ref(configuration.atp_maintenance_sbos), rxn_strings),
        )
        add_interface(:demand, sbod.(Ref(configuration.demand_sbos), rxn_strings))
    elseif interface == :identifier_prefixes
        prefixed(ps, s) = any(p -> startswith(s, p), ps)
        add_interface(
            :exchanges,
            prefixed.(Ref(configuration.exchange_id_prefixes), rxn_strings),
        )
        add_interface(
            :biomass,
            prefixed.(Ref(configuration.biomass_id_prefixes), rxn_strings),
        )
        add_interface(
            :atp_maintenance,
            in.(rxn_strings, Ref(configuration.atp_maintenance_ids)),
        )
    elseif interface == :boundary
        add_interface(
            :boundary,
            [(all(col .<= 0) | all(col .>= 0)) for col in eachcol(stoi)],
        )
    else
        interface == nothing || throw(DomainError(interface, "unknown interface specifier"))
    end

    return constraints
end

export flux_balance_constraints

"""
$(TYPEDSIGNATURES)

Compute an optimal objective-optimizing solution of the given `model`.

Most arguments are forwarded to [`optimized_values`](@ref).

Returns a tree with the optimization solution of the same shape as
given by [`flux_balance_constraints`](@ref).
"""
flux_balance_analysis(model::A.AbstractFBCModel; kwargs...) = frontend_optimized_values(
    flux_balance_constraints,
    model;
    objective = x -> x.objective.value,
    kwargs...,
)

export flux_balance_analysis
