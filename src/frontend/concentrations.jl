
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

Build log-concentration-stoichiometry constraints for the `model`, as used e.g.
by [`max_min_driving_force_analysis`](@ref).

The output constraint tree contains a log-concentration variable for each
`metabolite` from `model`, in subtree `log_concentrations`. The total reactant
log-concentrations for each `reaction` are constrained in subtree
`log_concentration_stoichiometry`. By default, all reactions and metabolites in
`model` are included.

A concentration bound is given by parameter function `concentration_bound` for
each metabolite ID (the string ID is the single argument of the function); by
default the function returns `nothing` and no bounds are installed. The same is
used for reactions with `reaction_concentration_bound`.
"""
function log_concentration_constraints(
    model::A.AbstractFBCModel;
    reactions = A.reactions(model),
    metabolites = A.metabolites(model),
    metabolite_concentration_bound = _ -> nothing,
    reaction_concentration_bound = _ -> nothing,
)
    rxns = String.(collect(reactions))
    mets = String.(collect(metabolites))
    metset = Set{String}(metabolites)

    vars = C.variables(keys = Symbol.(mets), bounds = metabolite_concentration_bound.(mets))

    stoi = C.ConstraintTree(
        Symbol(rxn) => C.Constraint(
            value = sum(
                (
                    coeff * vars[Symbol(met)].value for
                    (met, coeff) in A.reaction_stoichiometry(model, rxn) if (met in metset)
                );
                init = zero(C.LinearValue),
            ),
            bound = reaction_concentration_bound(rxn),
        ) for rxn in rxns
    )

    return :log_concentrations^vars * :log_concentration_stoichiometry^stoi
end

export log_concentration_constraints
