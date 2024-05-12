
# This file tests stuff that is too boring or mechanical to be included in the
# documentation. If you want to add tests here, first consider actually
# documenting the functionality in `docs/src/examples/` instead.

@testset "Switch bound" begin
    x = Switch(5, 10)
    y = -(((0.5 * (((1 - x) + 1) * 4)) / 2) - 2)
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
        optimizer = GLPK.Optimizer,
        settings = [set_objective_sense(Minimal), set_optimizer(GLPK.Optimizer)],
    ).x == 0.0
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
