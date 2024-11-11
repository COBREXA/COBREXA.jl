
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

Clamp all negative values in the bound to zero, leaving only the "positive
contribution" to the overall value of the constraint. Used in
[`unsigned_positive_contribution_variables`](@ref) and
[`unsigned_negative_contribution_variables`](@ref) to allocate unidirectional
variables.
"""
positive_bound_contribution(b::C.EqualTo) = b.equal_to >= 0 ? b : C.EqualTo(0.0)

positive_bound_contribution(b::C.Between) =
    b.lower >= 0 && b.upper >= 0 ? b :
    b.lower <= 0 && b.upper <= 0 ? C.EqualTo(0) :
    C.Between(max(0, b.lower), max(0, b.upper))

positive_bound_contribution(b::Switch) =
    let upper_bound = max(b.a, b.b)
        upper_bound > 0 ? C.Between(0.0, upper_bound) : C.EqualTo(0.0)
    end

"""
$(TYPEDSIGNATURES)

A constraint tree of variables with positive contributions to the values in `cs`.
"""
unsigned_positive_contribution_variables(cs::C.ConstraintTree) =
    C.variables_for(cs) do c
        positive_bound_contribution(something(c.bound, C.Between(-Inf, Inf)))
    end

export unsigned_positive_contribution_variables

"""
$(TYPEDSIGNATURES)

A constraint tree of variables with negative contributions to the values in `cs`.
"""
unsigned_negative_contribution_variables(cs::C.ConstraintTree) =
    C.variables_for(cs) do c
        positive_bound_contribution(-something(c.bound, C.Between(-Inf, Inf)))
    end


export unsigned_negative_contribution_variables

"""
$(TYPEDSIGNATURES)

A constraint tree that connects negative unsigned variable contributions to
signed ones, while acting as positive contributions.
"""
unsigned_positive_contribution_constraints(
    cs::C.ConstraintTree,
    negative::C.ConstraintTree,
) =
    C.zip(cs, negative) do c, n
        C.Constraint(
            c.value + n.value,
            positive_bound_contribution(something(c.bound, C.Between(-Inf, Inf))),
        )
    end

export unsigned_positive_contribution_constraints

"""
$(TYPEDSIGNATURES)

A constraint tree that connects positive unsigned variable contributions to
signed ones, while acting as negative contributions.
"""
unsigned_negative_contribution_constraints(
    cs::C.ConstraintTree,
    positive::C.ConstraintTree,
) =
    C.zip(cs, positive) do c, p
        C.Constraint(
            p.value - c.value,
            positive_bound_contribution(-something(c.bound, C.Between(-Inf, Inf))),
        )
    end

export unsigned_negative_contribution_constraints

"""
$(TYPEDSIGNATURES)

Shortcut for making a pair of named variable groups created by
[`unsigned_positive_contribution_variables`](@ref) and
[`unsigned_negative_contribution_variables`](@ref), in subtrees named by
`positive` and `negative`.

Use [`sign_split_constraints`](@ref) to bind the new variables to existing
values.
"""
sign_split_variables(
    constraints::C.ConstraintTree;
    positive::Symbol = :positive,
    negative::Symbol = :negative,
)::C.ConstraintTree =
    positive^unsigned_positive_contribution_variables(constraints) +
    negative^unsigned_negative_contribution_variables(constraints)

export sign_split_variables

"""
$(TYPEDSIGNATURES)

A constraint tree that bound the values present in `signed` to be sums of pairs
of `positive` and `negative` contributions to the individual values.

Keys in the result are the same as the keys of `signed` constraints.

Typically, this can be used to create "unidirectional" fluxes
together with [`unsigned_negative_contribution_variables`](@ref) and
[`unsigned_positive_contribution_variables`](@ref).

Use [`sign_split_variables`](@ref) to allocate the variables easily.
"""
sign_split_constraints(;
    positive::C.ConstraintTree,
    negative::C.ConstraintTree,
    signed::C.ConstraintTree,
) =
    C.zip(positive, negative, signed, C.Constraint) do p, n, s
        equal_value_constraint(s.value + n.value, p.value)
    end

export sign_split_constraints
