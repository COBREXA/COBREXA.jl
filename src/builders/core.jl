
"""
$(TYPEDSIGNATURES)

A constraint tree that models the content of the given instance of
`AbstractFBCModel`.
"""
function fbc_model_constraints(model::A.AbstractFBCModel)
    rxns = Symbol.(A.reactions(model))
    mets = Symbol.(A.metabolites(model))
    lbs, ubs = A.bounds(model)
    stoi = A.stoichiometry(model)
    bals = A.balance(model)
    obj = A.objective(model)

    #TODO: is sparse() required below?
    return C.ConstraintTree(
        :fluxes => C.variables(keys = rxns, bounds = zip(lbs, ubs)),
        :balances => C.ConstraintTree(
            m => C.Constraint(value = C.LinearValue(sparse(row)), bound = b) for
            (m, row, b) in zip(mets, eachrow(stoi), bals)
        ),
        :objective => C.Constraint(value = C.LinearValue(sparse(obj))),
    )
end

export fbc_model_constraints

"""
$(TYPEDSIGNATURES)

Shortcut for allocation non-negative ("unsigned") variables. The argument
`keys` is forwarded to `ConstraintTrees.variables` as `keys`.
"""
unsigned_variables(; keys) = C.variables(; keys, bounds = Ref((0.0, Inf)))

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
        value = s + (haskey(negative, k) ? C.value(negative[k]) : zero(C.Value)) -
                (haskey(positive, k) ? C.value(positive[k]) : zero(C.Value)),
        bound = 0.0,
    ) for (k, s) in C.elems(signed)
)
#TODO the example above might as well go to docs

export sign_split_constraints
