
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

# This file tests stuff that is too boring or mechanical to be included in the
# documentation. If you want to add tests here, first consider actually
# documenting the functionality in `docs/src/examples/` instead.

@testset "Switch bound" begin
    x = Switch(5, 10)
    y = -(((1 + 0.5 * (((1 - x) + 1) * 4)) / 2) - 2.5)
    @test isapprox(x.a, y.a)
    @test isapprox(x.b, y.b)
end

@testset "Sampling fails where it should" begin
    @test_throws DomainError sample_constraints(
        sample_chain_achr,
        :cons^C.Constraint(C.LinearValue([1], [1]), Switch(0, 1)),
        start_variables = fill(1.0, 1, 1),
        seed = UInt(1),
        n_chains = 1,
    )
end

@testset "Miscellaneous settings" begin
    c = :x^C.variable(bound = (0, 1))
    @test optimized_values(
        c,
        objective = c.x.value,
        optimizer = HiGHS.Optimizer,
        settings = [set_objective_sense(Minimal), set_optimizer(HiGHS.Optimizer)],
    ).x == 0.0
end

@testset "Builders" begin
    c = C.variables(keys = [:x, :y], bounds = C.Between(0, 1))
    a = less_or_equal_constraint(c.x, c.y)
    b = -greater_or_equal_constraint(c.y, c.x)
    @test a.value.idxs == b.value.idxs && a.value.weights == b.value.weights
    @test a.bound.lower == b.bound.lower && a.bound.upper == b.bound.upper
end

@testset "Bounds" begin
    x = COBREXA.positive_bound_contribution(Switch(1.0, 2.0))
    @test typeof(x) == C.Between
    @test x.lower == 0.0
    @test x.upper == 2.0
    x = COBREXA.positive_bound_contribution(C.EqualTo(-1.0))
    @test typeof(x) == C.EqualTo
    @test x.equal_to == 0.0
    x = scale_bounds(:a^C.variable(bound = (-5.0, 5.0)), 2)
    @test typeof(x.a.bound) == C.Between
    @test x.a.bound.lower == -10.0
    @test x.a.bound.upper == 10.0
end

@testset "Failing parsimonious objectives" begin
    c = :x^C.variable(bound = (0, 1))
    r = parsimonious_optimized_values(
        c,
        objective = c.x.value,
        parsimonious_objective = c.x.value,
        objective_value = 1.015,
        tolerances = [absolute_tolerance_bound(0.01 * i) for i = 1:5],
        optimizer = HiGHS.Optimizer,
    )
    @test r.x > 0.99
    r = parsimonious_optimized_values(
        c,
        objective = c.x.value,
        parsimonious_objective = c.x.value,
        objective_value = 1.0001,
        tolerances = [absolute_tolerance_bound(0.00001)],
        optimizer = HiGHS.Optimizer,
    )
    @test isnothing(r)
end

@testset "Many ways to specify an interface" begin
    x = C.variables(keys = [:x, :y], bounds = (0, 1))
    x *= :interface^x
    c = interface_constraints([
        :a => x,
        :b => (x, 2),
        :c => (x, :interface),
        :d => (x, :interface, 3),
    ])

    @test issetequal(keys(c.interface), [:x, :y])
    @test c.interface_balance.x.value.idxs == [1, 3, 5, 7, 9]
    @test c.interface_balance.x.value.weights == [1.0, 2.0, 1.0, 3.0, -1.0]

end

@testset "Uncommon bounds in gapfilling" begin
    vs = C.variables(keys = [:a, :b], bounds = [C.EqualTo(123), nothing])
    stoi = C.variables(keys = [:c, :d], bounds = C.EqualTo(0))
    x = gap_filling_constraints(
        system = vs,
        stoichiometry = stoi,
        universal_fluxes = vs,
        universal_stoichiometry = stoi,
    )
    @test x.universal_flux_bounds.a.bound isa C.EqualTo
    @test isempty(x.universal_flux_bounds.b)
    vs = C.variables(keys = [:a, :b], bounds = [C.EqualTo(123), Switch(1, 2)])
    @test_throws DomainError gap_filling_constraints(
        system = vs,
        stoichiometry = stoi,
        universal_fluxes = vs,
        universal_stoichiometry = stoi,
    )
end

@testset "Enzyme capacity expansion compat & corner cases" begin
    x = C.EqualTo(123.0)
    all = [:ident]
    y = expand_enzyme_capacity(x, all)
    @test length(y) == 1
    (_, (ks, v)) = y[1]
    # these things should not be touched, thus triple =
    @test ks === ident
    @test v === x

    x = expand_enzyme_capacity([(:test, [:ident], 123)], [:defa, :ults])
    @test length(x) == 1
    (i, (ks, v)) = x
    @test i == :test
    @test ks == [:ident]
    @test v isa C.Between
    @test v.lower == 0.0
    @test v.upper == 123.0
end
