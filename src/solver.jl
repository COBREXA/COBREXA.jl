
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

Representation of a "binary switch" bound for `ConstraintTree`s. The value is
constrained to be either the value of field `a` or of field `b`; both fields
are `Float64`s. Upon translation to JuMP, the switches create an extra boolean
variable, and the value is constrained to equal `a + boolean_var * (b-a)`.

Switches can be offset by adding real numbers, negated, and multiplied and
divided by scalar constraints. For optimizing some special cases, multiplying
by exact zero returns an equality bound to zero.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Switch <: C.Bound
    "One choice"
    a::Float64

    "The other choice"
    b::Float64
end

export Switch

Base.:-(x::Switch) = Switch(-s.a, -s.b)
Base.:+(x::Real, s::Switch) = b + a
Base.:+(s::Switch, x::Real) = Switch(s.a + x, s.b + x)
Base.:*(x::Real, s::Switch) = b * a
Base.:*(s::Switch, x::Real) = x == 0 ? C.EqualTo(0) : Switch(s.a * x, s.b * x)
Base.:/(s::Switch, x::Real) = Switch(s.a / x, s.b / x)

"""
$(TYPEDSIGNATURES)

Very efficiently substitute a ConstraintTrees' `LinearValue` into a JuMP
expression of type `AffExpr`.
"""
function substitute_jump(val::C.LinearValue, vars)
    e = J.AffExpr() # unfortunately @expression(model, 0) is not type stable and gives an Int
    for (i, w) in zip(val.idxs, val.weights)
        if i == 0
            J.add_to_expression!(e, w)
        else
            J.add_to_expression!(e, w, vars[i])
        end
    end
    e
end

"""
$(TYPEDSIGNATURES)

Very efficiently substitute a ConstraintTrees' `QuadraticValue` into a JuMP
expression of type `QuadExpr`.
"""
function substitute_jump(val::C.QuadraticValue, vars)
    e = J.QuadExpr()
    for ((i, j), w) in zip(val.idxs, val.weights)
        if i == 0 && j == 0
            J.add_to_expression!(e, w)
        elseif i == 0 # the symmetric case is prohibited
            J.add_to_expression!(e, w, vars[j])
        else
            J.add_to_expression!(e, w, vars[i], vars[j])
        end
    end
    e
end

"""
$(TYPEDSIGNATURES)

Add an equality constraint to a JuMP model.
"""
function constraint_jump!(model, expr, b::C.EqualTo)
    J.@constraint(model, expr == b.equal_to)
end

"""
$(TYPEDSIGNATURES)

Add an interval constraint to a JuMP model.
"""
function constraint_jump!(model, expr, b::C.Between)
    isinf(b.lower) || J.@constraint(model, expr >= b.lower)
    isinf(b.upper) || J.@constraint(model, expr <= b.upper)
end

"""
$(TYPEDSIGNATURES)

Add a [`Switch`](@ref) constraint to a JuMP model.
"""
function constraint_jump!(model, expr, b::Switch)
    boolean = J.@variable(model, binary = true)
    J.@constraint(model, expr == b.a + boolean * (b.b - b.a))
end

"""
$(TYPEDSIGNATURES)

Add an empty constraint to a JuMP model (this does not do anything).
"""
constraint_jump!(_, _, _::Nothing) = nothing

"""
$(TYPEDSIGNATURES)

Construct a JuMP `Model` that describes the precise constraint system into the
JuMP `Model` created for solving in `optimizer`, with a given optional
`objective` and optimization `sense` chosen from [`Maximal`](@ref),
[`Minimal`](@ref) and [`Feasible`](@ref).

All types of values in the constraint tree must have an overload for
[`substitute_jump`](@ref).
"""
function optimization_model(
    cs::C.ConstraintTreeElem;
    objective::Union{Nothing,C.Value} = nothing,
    optimizer,
    sense = Maximal,
)
    model = J.Model(optimizer)

    J.@variable(model, x[1:C.var_count(cs)])
    isnothing(objective) || J.@objective(model, sense, substitute_jump(objective, x))

    C.traverse(cs) do c
        isnothing(c.bound) || constraint_jump!(model, substitute_jump(c.value, x), c.bound)
    end

    return model
end

export optimization_model

"""
$(TYPEDSIGNATURES)

`true` if `opt_model` solved successfully (solution is optimal or
locally optimal). `false` if any other termination status is reached.
"""
is_solved(opt_model::J.Model) =
    J.termination_status(opt_model) in [J.MOI.OPTIMAL, J.MOI.LOCALLY_SOLVED]

export is_solved

"""
    Minimal

Objective sense for finding the minimal value of the objective.

Same as `JuMP.MIN_SENSE`.
"""
const Minimal = J.MIN_SENSE
export Minimal

"""
    Maximal

Objective sense for finding the maximal value of the objective.

Same as `JuMP.MAX_SENSE`.
"""
const Maximal = J.MAX_SENSE
export Maximal

"""
    Maximal

Objective sense for finding the any feasible value of the objective.

Same as `JuMP.FEASIBILITY_SENSE`.
"""
const Feasible = J.FEASIBILITY_SENSE
export Feasible

"""
$(TYPEDSIGNATURES)

Like [`optimized_values`](@ref), but works directly with a given JuMP
model `om` without applying any settings or creating the optimization model.

To run the process manually, you can use [`optimization_model`](@ref) to
convert the constraints into a suitable JuMP optimization model.
"""
function optimized_model(om; output::C.ConstraintTreeElem)
    J.optimize!(om)
    is_solved(om) ? C.substitute_values(output, J.value.(om[:x])) : nothing
end

export optimized_model

"""
$(TYPEDSIGNATURES)

Like [`optimized_model`](@ref) but only returns the objective value (or
`nothing` if the model is not solved).
"""
function optimized_objective(om)
    J.optimize!(om)
    is_solved(om) ? J.objective_value(om) : nothing
end
