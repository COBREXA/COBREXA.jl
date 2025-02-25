
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

Join multiple constraint tree modules with interfaces into a bigger module with
an interface.

Modules are like usual constraint trees, but contain an explicitly declared
`interface` part, marked properly in arguments using e.g. a tuple (the
parameters should form a dictionary constructor that would generally look such
as `:module_name => (module, module.interface)`; the second tuple member may
also be specified just by name as e.g. `:interface`, or omitted while relying
on `default_interface`).

Interface parts get merged and constrained to create a new interface; networks
are intact with disjoint variable sets.

Compatible modules with ready-made interfaces may be created e.g. by
[`flux_balance_constraints`](@ref).

`ignore` may be used to selectively ignore parts of interfaces given the
"module name" identifier and constraint path in the interface (these form 2
parameters passed to `ignore`). Similarly, `bound` may be used to specify
bounds for the new interface, if required.

`output_modules` is a function to transform the contraint tree of the assembled
modules before adding the connection structure; this can be used to e.g. put
the members into a single subtree with `output_modules = x -> :subtree_name^x`.
"""
function interface_constraints(
    ps::Pair...;
    default_interface = :interface,
    out_interface = default_interface,
    out_balance = Symbol(out_interface, :_balance),
    output_modules = identity,
    ignore = (_, _) -> false,
    bound = _ -> nothing,
)

    prep(id::String, x) = prep(Symbol(id), x)
    prep(id::Symbol, mod::C.ConstraintTree) = prep(id, (mod, default_interface))
    prep(id::Symbol, (mod, multiplier)::Tuple{C.ConstraintTree,<:Real}) =
        prep(id, (mod, default_interface, multiplier))
    prep(id::Symbol, (mod, interface)::Tuple{C.ConstraintTree,Symbol}) =
        prep(id, (mod, mod[interface]))
    prep(id::Symbol, (mod, interface, multiplier)::Tuple{C.ConstraintTree,Symbol,<:Real}) =
        prep(id, (mod, mod[interface], multiplier))
    prep(
        id::Symbol,
        (mod, interface, multiplier)::Tuple{C.ConstraintTree,C.ConstraintTreeElem,<:Real},
    ) = prep(id, (mod, scale_constraints(interface, multiplier)))
    prep(id::Symbol, (mod, interface)::Tuple{C.ConstraintTreeElem,C.ConstraintTreeElem}) =
        (id^(:network^mod * :interface^interface))
    prep_pair((a, b)) = prep(a, b)

    # first, collect everything into one huge network
    # (while also renumbering the interfaces)
    modules = sum(prep_pair.(ps); init = C.ConstraintTree())

    # TODO maybe split the interface-preparing function into a separate
    # function, so that others can easily make their own interface connectors
    # with compatible interface.

    # fold a union of all non-ignored interface keys
    interface_sum = foldl(modules, init = C.ConstraintTree()) do accs, (id, ms)
        C.imerge(accs, ms.interface) do path, acc, m
            ignore(id, path) ? missing :
            ismissing(m) ? acc :
            ismissing(acc) ? C.Constraint(value = m.value) :
            C.Constraint(value = acc.value + m.value)
        end
    end

    # extract the plain networks and add variables for the new interfaces
    constraints =
        output_modules(C.ConstraintTree(id => (m.network) for (id, m) in modules)) +
        out_interface^C.variables_ifor((path, _) -> bound(path), interface_sum)

    # join everything with the interface balance and return
    constraints * out_balance^C.zip(interface_sum, constraints[out_interface]) do sum, out
        C.Constraint(sum.value - out.value, 0)
    end
end

"""
$(TYPEDSIGNATURES)

Overload of [`interface_constraints`](@ref) for general key-value containers.
"""
interface_constraints(kv; kwargs...) = interface_constraints(kv...; kwargs...)

export interface_constraints

"""
$(TYPEDSIGNATURES)

Glue an interface made of additional variables into a constraint system that
lacks variables to support one. Returns a tuple of a "module" and its
"interface", usable in [`interface_constraints`](@ref).

First, `interface` variables are renumbered to not collide with `constraints`.
The values in `interface` are then multiplied by `multiplier` (by default 1)
and injected into `constraints` that have the same "path" in the constraint
tree. Attempts to inject into missing `constraints` and mismatches between leaf
constraints and subtrees are ignored. The output consists of the extended
`constraints` and the renumbered `interface` (in particular, numbering of
the original variables in `constraints` is retained). Constraints in both trees
retain their bounds.

In metabolic modeling terms, this adds interfaces that contribute to any
constrained quantity in a model. For example, one may conveniently create new
exchange reactions to contribute to an existing metabolite balance, which can
in turn be used as an interface in a more complex model.
"""
function inject_interface(
    constraints::C.ConstraintTree,
    interface::C.ConstraintTree;
    multiplier = 1,
)
    go(cs::C.ConstraintTree, iface::C.ConstraintTree) =
        C.ConstraintTree(k => haskey(iface, k) ? go(v, iface[k]) : v for (k, v) in cs)
    go(cs::C.Constraint, iface::C.Constraint) =
        C.Constraint(value = cs.value + multiplier * iface.value, bound = cs.bound)
    go(cs, _) = cs # any non-matching interface is simply ignored

    interface_ = C.incr_var_idxs(interface, C.var_count(constraints))
    return (go(constraints, interface_), interface_)
end

export inject_interface
