
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

Make a function that returns absolute tolerance bounds, i.e. `value -
tolerance` and `value + tolerance` in a tuple, in the increasing order.
"""
absolute_tolerance_bound(tolerance) = x -> begin
    bound = (x - tolerance, x + tolerance)
    (minimum(bound), maximum(bound))
end

export absolute_tolerance_bound

"""
$(TYPEDSIGNATURES)

Make a function that returns relative tolerance bounds, i.e. `value /
tolerance` and `value * tolerance` in a tuple, in the increasing order.
"""
relative_tolerance_bound(tolerance) = x -> begin
    bound = (x * tolerance, x / tolerance)
    (minimum(bound), maximum(bound))
end

export relative_tolerance_bound

"""
$(TYPEDSIGNATURES)

Make a copy of a constraint tree with all bounds removed. This is helpful when
creating large trees only for for value representation purposes, which should
not directly constraint anything (and thus should not put additional stress on
the constraint solver).
"""
remove_bounds(cs::C.ConstraintTree) =
    C.map(cs) do c
        C.Constraint(c.value, nothing)
    end

export remove_bounds
