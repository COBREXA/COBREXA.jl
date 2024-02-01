
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

Create max-min driving force constraint system from `model`, using the supplied
reaction standard Gibbs free energies in `reaction_standard_gibbs_free_energies`.

The method is described by: *Noor, et al., "Pathway thermodynamics highlights
kinetic obstacles in central metabolism.", PLoS computational biology, 2014.*

`reference_flux` sets the directions of each reaction in `model`. The scale of
the values is not important, only the direction is examined (w.r.t.
`reference_flux_atol` tolerance). Ideally, the supplied `reference_flux` should
be completely free of internal cycles, which enables the thermodynamic
consistency. To get the cycle-free flux, you can use
[`loopless_flux_balance_analysis`](@ref) (computationally demanding, but gives
thermodynamically consistent solutions),
[`parsimonious_flux_balance_analysis`](@ref) or
[`linear_parsimonious_flux_balance_analysis`](@ref) (which is computationally
simple, but the consistency is not guaranteed).

Internally, [`log_concentration_constraints`](@ref) is used to lay out the
base structure of the problem.

Following arguments are set optionally:
- `water_metabolites`, `proton_metabolites` and `ignored_metabolites` allow to
  completely ignore constraints on a part of the metabolite set, which is
  explicitly recommended especially for water and protons (since the analyses
  generally assume aqueous environment of constant pH)
- `constant_concentrations` can be used to fix the concentrations of the
  metabolites
- `concentration_lower_bound` and `concentration_upper_bound` set the default
  concentration bounds for all other metabolites
- `concentration ratios` is a dictionary that assigns a tuple of
  metabolite-metabolite-concentration ratio constraint names; e.g. ATP/ADP
  ratio can be fixed to five-times-more-ATP by setting `concentration_ratios =
  Dict("adenosin_ratio" => ("atp", "adp", 5.0))`
- `T` and `R` default to the usual required thermodynamic constraints in the
  usual units (K and kJ/K/mol, respectively). These multiply the
  log-concentrations to obtain the actual Gibbs energies and thus driving
  forces.
"""
function max_min_driving_force_constraints(
    model::A.AbstractFBCModel;
    reaction_standard_gibbs_free_energies::Dict{String,Float64},
    reference_flux = Dict{String,Float64}(),
    concentration_ratios = Dict{String,Tuple{String,String,Float64}}(),
    constant_concentrations = Dict{String,Float64}(),
    ignored_metabolites = [],
    proton_metabolites = [],
    water_metabolites = [],
    concentration_lower_bound = 1e-9, # Molar
    concentration_upper_bound = 1e-1, # Molar
    T = 298.15, # Kelvin
    R = 8.31446261815324e-3, # kJ/K/mol
    reference_flux_atol = 1e-6,
)

    #=
    For MMDF, one is almost always interested in running the analysis on a
    subset of the full model, because thermodynamic data is usually not
    available for all the reactions/metabolites. Missing data may cause subtle
    problems, and is usually best to restrict your model to reactions where
    thermodynamic data exists.  
    =#

    model_reactions =
        isempty(reference_flux) ? A.reactions(model) :
        collect(k for (k, v) in reference_flux if abs(v) >= reference_flux_atol)

    model_metabolites =
        isempty(reference_flux) ? A.metabolites(model) :
        unique(
            mid for rid in model_reactions for
            mid in keys(A.reaction_stoichiometry(model, string(rid)))
        )

    # First let's just check if all the identifiers are okay because we use
    # quite a lot of these; the rest of the function may be a bit cleaner with
    # this checked properly.


    all(in(keys(reaction_standard_gibbs_free_energies)), model_reactions) || throw(
        DomainError(
            reaction_standard_gibbs_free_energies,
            "missing reaction in reaction_standard_gibbs_free_energies",
        ),
    )
    all(in(model_reactions), keys(reference_flux)) || throw(
        DomainError(
            reaction_standard_gibbs_free_energies,
            "unknown reactions referenced by reference_flux",
        ),
    )
    all(in(model_metabolites), keys(constant_concentrations)) || throw(
        DomainError(
            constant_concentrations,
            "unknown metabolites referenced by constant_concentrations",
        ),
    )
    all(
        in(model_metabolites),
        (m for (_, (x, y, _)) in concentration_ratios for m in (x, y)),
    ) || throw(
        DomainError(
            concentration_ratios,
            "unknown metabolites referenced by concentration_ratios",
        ),
    )
    all(in(model_metabolites), proton_metabolites) || throw(
        DomainError(
            concentration_ratios,
            "unknown metabolites referenced by proton_metabolites",
        ),
    )
    all(in(model_metabolites), water_metabolites) || throw(
        DomainError(
            concentration_ratios,
            "unknown metabolites referenced by water_metabolites",
        ),
    )
    all(in(model_metabolites), ignored_metabolites) || throw(
        DomainError(
            concentration_ratios,
            "unknown metabolites referenced by ignored_metabolites",
        ),
    )

    # that was a lot of checking.

    default_concentration_bound =
        C.Between(log(concentration_lower_bound), log(concentration_upper_bound))

    no_concentration_metabolites = union(
        Set(Symbol.(water_metabolites)),
        Set(Symbol.(proton_metabolites)),
        Set(Symbol.(ignored_metabolites)),
    )

    # TODO it might be quite viable to only create a system with the
    # metabolites and reactions actually required for the MMDF.

    constraints =
        log_concentration_constraints(
            model,
            subset_reactions = rxn -> rxn in model_reactions,
            subset_metabolites = met -> met in model_metabolites,
            concentration_bound = met -> if Symbol(met) in no_concentration_metabolites
                C.EqualTo(0)
            else
                if haskey(constant_concentrations, met)
                    C.EqualTo(log(constant_concentrations[met]))
                else
                    default_concentration_bound
                end
            end,
        ) + :min_driving_force^C.variable()

    driving_forces = C.ConstraintTree(
        let r = Symbol(rid),
            rf = reference_flux[rid],
            dGr0 = reaction_standard_gibbs_free_energies[rid],
            df = dGr0 + R * T * constraints.reaction_stoichiometry[r].value

            r => C.Constraint(rf >= 0 ? -df : df, C.Between(0, Inf))

        end for rid in model_reactions
    )

    constraints *
    :driving_forces^driving_forces *
    :min_driving_force_thresholds^C.map(driving_forces) do c
        less_or_equal_constraint(constraints.min_driving_force, c)
    end *
    :concentration_ratio_constraints^C.ConstraintTree(
        Symbol(cid) => difference_constraint(
            constraints.log_concentrations[Symbol(m1)],
            constraints.log_concentrations[Symbol(m2)],
            log(ratio),
        ) for (cid, (m1, m2, ratio)) in concentration_ratios
    )
end

export max_min_driving_force_constraints

"""
$(TYPEDSIGNATURES)

Perform the max-min driving force analysis on the `model`, returning the solved
constraint system.

Arguments are forwarded to [`max_min_driving_force_constraints`](@ref) (see the
documentation for the description of the constraint system); solver
configuration arguments are forwarded to [`optimized_values`](@ref).
"""
max_min_driving_force_analysis(model::A.AbstractFBCModel; kwargs...) =
    frontend_optimized_values(
        max_min_driving_force_constraints,
        model;
        objective = x -> x.min_driving_force.value,
        kwargs...,
    )

export max_min_driving_force_analysis
