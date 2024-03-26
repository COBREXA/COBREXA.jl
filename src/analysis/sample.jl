
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
function sampler_chain_achr(;
    sample::M;
    single_variable_lower_bounds::V,
    single_variable_upper_bounds::V,
    eq_constraints::M,
    eq_values::V,
    interval_constraints::M,
    interval_lower_bounds::V,
    interval_upper_bounds::V,
    eq_epsilon::F,
    interval_epsilon::F,
    collect_iterations::Vector{Int},
    generator::StableRNG.StableRNG,
) where {F<:Real,V<:AbstractVector{F},M<:AbstractMatrix{F}}

    (n, d) = size(sample)
    nec = size(eq_constraints, 1)
    nic = size(interval_constraints, 1)

    @assert length(single_variable_lower_bounds) == d
    @assert length(single_variable_upper_bounds) == d
    @assert size(eq_constraints, 2) == d
    @assert nec == length(eq_values)
    @assert size(interval_constraints, 2) == d
    @assert nic == length(interval_lower_bounds) == length(interval_upper_bounds)

    result = Matrix{F}(undef, d, n * length(collect_iterations))
    center = Vector{F}(undef, d)

    function update_range(range, pos, dir, lb, ub, tol)
        (rl, ru) = range
        dl = lb - pos
        du = ub - pos
        lower, upper =
            dir < -tol ? (du, dl) ./ dir : dir > tol ? (dl, du) ./ dir : (F(-Inf), F(Inf))
        return (max(rl, lower), min(ru, upper))
    end

    iter = 0

    for (iter_idx, iter_target) in enumerate(collect_iterations)
        while iter < iter_target
            iter += 1

            center .= sum(sample, dims = 2) ./ n

            for i = 1:n
                dir = center .- sample[:, i]

                # iteratively collect the maximum and minimum possible multiple
                # of `dir` added to the current point
                run_range = (F(-Inf), F(Inf))
                for j = 1:d
                    run_range = update_range(
                        run_range,
                        sample[j, i],
                        dir[j],
                        single_variable_lower_bounds[j],
                        single_variable_upper_bounds[j],
                        interval_epsilon,
                    )
                end

                # do the same for coupling constraints (equality first)
                dc = eq_constraints * dir
                pc = eq_constraints * sample[:, i]
                for j = 1:nec
                    run_range = update_range(
                        run_range,
                        pc[j],
                        dc[j],
                        eq_values[j],
                        eq_values[j],
                        eq_epsilon,
                    )
                end

                # ... and intervals
                dc = interval_constraints * dir
                pc = interval_constraints * sample[:, i]
                for j = 1:nic
                    run_range = update_range(
                        run_range,
                        pc[j],
                        dc[j],
                        interval_lower_bounds[j],
                        interval_upper_bounds[j],
                        interval_epsilon,
                    )
                end

                # generate a point in the viable run range and update it
                lambda = run_range[1] + rand(rng) * (run_range[2] - run_range[1])
                isfinite(lambda) || continue # bail out to avoid divergence
                sample[:, i] .= points[:, i] .+ lambda .* dir
            end
        end
        result[:, n*(iter_idx-1)+1:iter_idx*n] .= sample
    end

    return sample
end

"""
$(TYPEDSIGNATURES)

TODO
"""
function sample_constraint_variables(
    sampler::Function,
    constraints::C.ConstraintTree,
    start_variables::Matrix{Float64},
)

    return sample
end

"""
$(TYPEDSIGNATURES)

TODO
"""
function sample_constraints(
    sampler::Function,
    constraints::C.ConstraintTree,
    output::C.ConstraintTree,
    aggregate = collect,
    sample::Matrix{Float64},
)

end
