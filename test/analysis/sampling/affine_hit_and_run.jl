@testset "Sampling Tests" begin

    model = load_model(StandardModel, model_paths["e_coli_core.json"])

    warmup = warmup_from_variability(model, Tulip.Optimizer; workers = W)

    samples = affine_hit_and_run(
        model,
        warmup,
        sample_iters = 10 * (1:3),
        workers = W,
        chains = length(W),
    )

    @test size(samples, 1) == size(warmup, 1)
    @test size(samples, 2) == size(warmup, 2) * 3 * length(W)

    lbs, ubs = bounds(model)
    @test all(samples .>= lbs)
    @test all(samples .<= ubs)
    @test all(stoichiometry(model) * samples .< TEST_TOLERANCE)
end
