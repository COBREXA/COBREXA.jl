
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
function sample_chain_achr(;
    sample::M;
    variable_lower_bounds::V,
    variable_upper_bounds::V,
    coupling::M,
    lower_bounds::V,
    upper_bounds::V,
    epsilon::F = config.sampler_tolerance,
    collect_iterations::Vector{Int},
    generator::StableRNG.StableRNG,
) where {F<:Real,V<:AbstractVector{F},M<:AbstractMatrix{F}}

    (n, d) = size(sample)
    nc = size(coupling, 1)

    @assert length(variable_lower_bounds) == d
    @assert length(variable_upper_bounds) == d
    @assert size(coupling, 2) == d
    @assert nc == length(lower_bounds) == length(upper_bounds)

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
                        variable_lower_bounds[j],
                        variable_upper_bounds[j],
                        epsilon,
                    )
                end

                # do the same for coupling constraints
                dc = coupling * dir
                pc = coupling * sample[:, i]
                for j = 1:nc
                    run_range = update_range(
                        run_range,
                        pc[j],
                        dc[j],
                        lower_bounds[j],
                        upper_bounds[j],
                        epsilon,
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
    constraints::C.ConstraintTree;
    start_variables::Matrix{Float64},
    seed::Int,
    workers = D.workers(),
    n_chains,
    kwargs...,
)
    # preallocate a bit of stuff
    d = size(start_variables, 2)
    @assert C.var_count(constraints) <= d
    lbs = fill(-Inf, d)
    ubs = fill(-Inf, d)
    cs = C.Constraint[]

    # gather constraints information
    tick(v::C.LinearValue, b::C.EqualTo) = tick(v, C.Between(b.equal_to, b.equal_to))
    tick(v::C.LinearValue, b::C.Between) =
        if length(v.idxs) == 1 && first(v.idxs) > 0
            idx = first(v.idxs)
            lbs[idx] = max(lbs[idx], b.lower)
            ubs[idx] = min(ubs[idx], b.lower)
        else
            push!(cs, C.Constraint(v, b))
        end
    tick(_, ::Nothing) = nothing
    tick(v, b) = throw(
        DomainError(
            (v, b),
            "sampling from unsupported constraint type: $(typeof(v)), $(typeof(b))",
        ),
    )

    C.traverse(constraints) do c
        tick(c.value, c.bound)
    end

    # convert constraints to a matrix
    coupling = vcat(transpose.(sparse(c.value) for c in cs)...)
    clbs = [c.bound.lower for c in cs]
    cubs = [c.bound.upper for c in cs]

    # run the actual sampling using screen()
    rng = StableRNG(seed)
    samples = screen(rand(rng, Int, n_chains)) do chain_seed
        sampler(
            start_variables;
            variable_lower_bounds = lbs,
            variable_upper_bounds = ubs,
            coupling,
            lower_bounds = clbs,
            upper_bounds = cubs,
            seed = chain_seed,
            kwargs...,
        )
    end
    return vcat(samples...)
end

"""
$(TYPEDSIGNATURES)

TODO
"""
function sample_constraints(
    sampler::Function,
    constraints::C.ConstraintTree;
    output::C.ConstraintTree = constraints,
    aggregate = collect,
    kwargs...,
)
    var_samples = [
        aggregate(col) for col in eachcol(
            sample_constraint_variables(sampler, constraints; start_variables, kwargs...),
        )
    ]

    C.map(output) do c
        C.substitute(c.value, var_samples)
    end
end
