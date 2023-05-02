
@testset "HDF5 model SBML model" begin
    model = load_model(MatrixModel, model_paths["e_coli_core.xml"])
    fn = "ecoli_test.h5.noextension"
    h5m = save_model(model, fn, extension = ".h5")
    @test h5m isa HDF5Model
    @test h5m.filename == fn
    @test h5m.h5 == nothing #the file should not be open by default

    h5 = load_model(fn, extension = ".h5")
    precache!(h5)
    @test !isnothing(h5.h5)

    # briefly test that the loading is okay
    @test variable_count(model) == variable_count(h5)
    @test n_metabolites(model) == n_metabolites(h5)
    @test issetequal(variable_ids(model), variable_ids(h5))
    @test issetequal(metabolites(model), metabolites(h5))
    @test issorted(metabolites(h5))
    @test issorted(variable_ids(h5))
    @test size(stoichiometry(model)) == size(stoichiometry(h5))
    @test isapprox(sum(stoichiometry(model)), sum(stoichiometry(h5)))
    rxnp = sortperm(variable_ids(model))
    @test bounds(model)[1][rxnp] == bounds(h5)[1]
    @test bounds(model)[2][rxnp] == bounds(h5)[2]
    @test objective(model)[rxnp] == objective(h5)
    @test all(iszero, balance(h5))

    close(h5)
    @test isnothing(h5.h5)
end
