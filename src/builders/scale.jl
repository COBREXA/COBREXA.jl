
# Copyright (c) 2021-2025, University of Luxembourg
# Copyright (c) 2021-2025, Heinrich-Heine University Duesseldorf
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

Linearly scale all bounds in a constraint tree by the `factor`. This actually
changes the model semantics, and may not work in surprising/improper ways with
some constraint systems, esp. the MILP and QP ones.

See also [`scale_constraints`](@ref).
"""
scale_bounds(tree::C.ConstraintTree, factor) =
    C.map(tree) do c
        isnothing(c.bound) ? c : C.Constraint(value = c.value, bound = factor * c.bound)
    end

export scale_bounds

"""
$(TYPEDSIGNATURES)

Linearly scale all constraints in a constraint tree by the `factor`.

See also [`scale_bounds`](@ref).
"""
scale_constraints(tree::C.ConstraintTree, factor) =
    C.map(tree) do c
        c * factor
    end

export scale_constraints

"""
$(TYPEDSIGNATURES)

Given the value `x` constrained by bound `b`, produce a constraint system where
the value is bounded by `b` scaled by `scale`. The result may contain multiple
constraints.

Additional overloads may be implemented to handle custom kinds of bounds.

This overload scales an equality constraint; producing a single constraint in
the form `x-bound*scale == 0`.
"""
value_scaled_bound_constraint(x::C.Value, b::C.EqualTo, scale::C.Value) =
    C.Constraint(x - b.equal_to * scale, 0)

"""
$(TYPEDSIGNATURES)

Overload of [`value_scaled_bound_constraint`](@ref) for interval bounds.

This produces up to 2 constraints:
- `lower` bound in the form `x + lower_bound * scale >= 0`
- `upper` bound in the form `x - upper_bound * scale <= 0`
"""
function value_scaled_bound_constraint(x::C.Value, b::C.Between, scale::C.Value)
    bounds = [
        b.lower > -Inf ? (:lower => C.Constraint(x - b.lower * scale, (0, Inf))) : nothing,
        b.upper < Inf ? (:upper => C.Constraint(x - b.upper * scale, (-Inf, 0))) : nothing,
    ]
    return C.ConstraintTree(b for b in bounds if !isnothing(b))
end

"""
$(TYPEDSIGNATURES)

No-op overload of [`value_scaled_bound_constraint`](@ref) for unbounded
constriants. Produces an empty constraint tree.
"""
value_scaled_bound_constraint(x::C.Value, b::Nothing, scale::C.Value) = C.ConstraintTree()

# TODO handle Switch bounds too

export value_scaled_bound_constraint

"""
$(TYPEDSIGNATURES)

Produce a constraint tree with all bounds scaled by the `scale`, as defined by
[`value_scaled_bound_constraint`](@ref).

Values that carry no bound (and generally all empty constraint trees in the
solution) are removed. To preserve the original values without the risk of
violating the scaled constraints, use [`remove_bounds`](@ref) on `x`.
"""
value_scaled_bound_constraints(x::C.ConstraintTree, scale::C.Value) = C.ConstraintTree(
    k => v for
    (k, v) in (k => value_scaled_bound_constraints(v, scale) for (k, v) in x) if
    !(v isa C.ConstraintTree) || !isempty(v)
)

"""
$(TYPEDSIGNATURES)

Convenience overload of [`value_scaled_bound_constraints`](@ref) for a single
constraint.
"""
value_scaled_bound_constraints(x::C.Constraint, scale::C.Value) =
    value_scaled_bound_constraint(x.value, x.bound, scale)

export value_scaled_bound_constraints
