
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

TODO
"""
loopless_constraints(;
    fluxes::C.ConstraintTree,
    loopless_direction_indicators::C.ConstraintTree,
    loopless_driving_forces::C.ConstraintTree,
    internal_reactions::Vector{Symbol},
    internal_nullspace::Matrix,
    flux_infinity_bound,
    driving_force_nonzero_bound,
    driving_force_infinity_bound,
) = C.ConstraintTree(
    :flux_direction_lower_bounds => C.ConstraintTree(
        r => C.Constraint(
            value = fluxes[r].value +
                    flux_infinity_bound * (1 - loopless_direction_indicators[r].value),
            bound = C.Between(0, Inf),
        ) for r in internal_reactions
    ),
    :flux_direction_upper_bounds => C.ConstraintTree(
        r => C.Constraint(
            value = fluxes[r].value +
                    flux_infinity_bound * loopless_direction_indicators[r].value,
            bound = C.Between(-Inf, 0),
        ) for r in internal_reactions
    ),
    :driving_force_lower_bounds => C.ConstraintTree(
        r => C.Constraint(
            value = loopless_driving_forces[r].value -
                    driving_force_nonzero_bound * loopless_direction_indicators[r].value +
                    driving_force_infinity_bound *
                    (1 - loopless_direction_indicators[r].value),
            bound = C.Between(0, Inf),
        ) for r in internal_reactions
    ),
    :driving_force_upper_bounds => C.ConstraintTree(
        r => C.Constraint(
            value = loopless_driving_forces[r].value +
                    driving_force_nonzero_bound *
                    (1 - loopless_direction_indicators[r].value) -
                    driving_force_infinity_bound * loopless_direction_indicators[r].value,
            bound = C.Between(-Inf, 0),
        ) for r in internal_reactions
    ),
    :loopless_nullspace => C.ConstraintTree(
        Symbol(:nullspace_base_, i) => C.Constraint(
            value = sum(
                coeff * loopless_driving_forces[r].value for
                (coeff, r) in zip(col, internal_reactions)
            ),
            bound = C.EqualTo(0),
        ) for (i, col) in enumerate(eachcol(internal_nullspace))
    ),
)

export loopless_constraints
