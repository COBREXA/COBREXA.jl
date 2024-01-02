
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
"""
function fbc_model_constraints(model::A.AbstractFBCModel)
    rxns = Symbol.(A.reactions(model))
    mets = Symbol.(A.metabolites(model))
    lbs, ubs = A.bounds(model)
    stoi = A.stoichiometry(model)
    bal = A.balance(model)
    obj = A.objective(model)

    #TODO: is sparse() required below?
    return C.ConstraintTree(
        :fluxes^C.variables(keys = rxns, bounds = zip(lbs, ubs)) *
        :flux_stoichiometry^C.ConstraintTree(
            met => C.Constraint(value = C.LinearValue(sparse(row)), bound = C.EqualTo(b))
            for (met, row, b) in zip(mets, eachrow(stoi), bal)
        ) *
        :objective^C.Constraint(C.LinearValue(sparse(obj))),
    )
end

export fbc_model_constraints

"""
$(TYPEDSIGNATURES)

Shortcut for allocation non-negative ("unsigned") variables. The argument
`keys` is forwarded to `ConstraintTrees.variables` as `keys`.
"""
unsigned_variables(; keys) = C.variables(; keys, bounds = C.Between(0.0, Inf))

export unsigned_variables

"""
$(TYPEDSIGNATURES)

A constraint tree that bound the values present in `signed` to be sums of pairs
of `positive` and `negative` contributions to the individual values.

Keys in the result are the same as the keys of `signed` constraints.

Typically, this can be used to create "unidirectional" fluxes
together with [`unsigned_variables`](@ref):
```
uvars = unsigned_variables(keys(myModel.fluxes))

myModel = myModel +
    :fluxes_forward^uvars +
    :fluxes_reverse^uvars

myModel *=
    :direction_sums^sign_split_constraints(
        positive = myModel.fluxes_forward,
        negative = myModel.fluxes_reverse,
        signed = myModel.fluxes,
    )
```
"""
sign_split_constraints(;
    positive::C.ConstraintTree,
    negative::C.ConstraintTree,
    signed::C.ConstraintTree,
) = C.ConstraintTree(
    k => C.Constraint(
        value = s.value +
                (haskey(negative, k) ? negative[k].value : zero(typeof(s.value))) -
                (haskey(positive, k) ? positive[k].value : zero(typeof(s.value))),
        bound = C.EqualTo(0.0),
    ) for (k, s) in signed
)
#TODO the example above might as well go to docs

export sign_split_constraints

function fluxes_in_direction(fluxes::C.ConstraintTree, direction = :forward)
    keys = Symbol[]
    for (id, flux) in fluxes
        if direction == :forward
            flux.bound.upper > 0 && push!(keys, id)
        else
            flux.bound.lower < 0 && push!(keys, id)
        end
    end
    C.variables(; keys, bounds = C.Between(0.0, Inf))
end

export fluxes_in_direction
