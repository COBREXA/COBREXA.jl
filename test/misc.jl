
# This file tests stuff that is too boring or mechanical to be included in the
# documentation. If you want to add tests here, first consider actually
# documenting the functionality in `docs/src/examples/` instead.

@testset "Switch bound" begin
    x = Switch(5, 10)
    y = -((((((-x) + 1) * 4) / 2) / 2) - 1)
    @test isapprox(x.a, y.a)
    @test isapprox(x.b, y.b)
end
