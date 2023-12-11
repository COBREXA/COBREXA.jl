"""
$(TYPEDSIGNATURES)

Allocate enzyme variables (gene products in the model) to a constraint tree.
"""
enzyme_variables(model::A.AbstractFBCModel) =
    C.variables(; keys = Symbol.(A.genes(model)), bounds = Ref((0.0, Inf)))

"""
$(TYPEDSIGNATURES)

Create isozyme variables for reactions. A single reaction may be catalyzed by
multiple enzymes (isozymes), and the total flux through a reaction is the sum
through of all these isozymes. These variables are linked to fluxes through
[`link_isozymes`](@ref).
"""
function isozyme_variables(
    reaction_id::String,
    reaction_isozymes::Dict{String,Dict{String,Isozyme}},
)
    C.variables(;
        keys = Symbol.(collect(keys(reaction_isozymes[reaction_id]))),
        bounds = Ref((0.0, Inf)),
    )
end

"""
$(TYPEDSIGNATURES)

Link isozymes to fluxes. All forward (backward) fluxes are the sum of all the
isozymes catalysing these fluxes.
"""
function link_isozymes(
    fluxes_directional::C.ConstraintTree,
    fluxes_isozymes::C.ConstraintTree,
)
    C.ConstraintTree(
        k => C.Constraint(
            value = s.value - sum(x.value for (_, x) in fluxes_isozymes[k]),
            bound = 0.0,
        ) for (k, s) in fluxes_directional if haskey(fluxes_isozymes, k)
    )
end

"""
$(TYPEDSIGNATURES)

Create the enzyme "mass balance" matrix. In essence, the stoichiometric
coefficient is subunit_stoichiometry / kcat for each directional, isozyme flux,
and it must be balanced by the enzyme variable supply.
"""
function enzyme_stoichiometry(
    enzymes::C.ConstraintTree,
    fluxes_isozymes_forward::C.ConstraintTree,
    fluxes_isozymes_backward::C.ConstraintTree,
    reaction_isozymes::Dict{String,Dict{String,Isozyme}},
)
    # map enzyme ids to reactions that use them (not all isozymes have to though)
    enzyme_rid_lookup = Dict{Symbol,Vector{Symbol}}()
    for (rid, isozymes) in reaction_isozymes
        for isozyme in values(isozymes)
            for gid in keys(isozyme.gene_product_stoichiometry)
                rids = get!(enzyme_rid_lookup, Symbol(gid), Symbol[])
                rid in rids || push!(rids, Symbol(rid))
            end
        end
    end

    C.ConstraintTree(
        gid => C.Constraint(
            value = enz.value + # supply
                    sum(
                        enzyme_balance(
                            gid,
                            rid,
                            fluxes_isozymes_forward,
                            reaction_isozymes,
                            :kcat_forward,
                        ) for rid in enzyme_rid_lookup[gid] if
                        haskey(fluxes_isozymes_forward, rid);
                        init = zero(typeof(enz.value)),
                    ) + # flux through positive isozymes
                    sum(
                        enzyme_balance(
                            gid,
                            rid,
                            fluxes_isozymes_backward,
                            reaction_isozymes,
                            :kcat_backward,
                        ) for rid in enzyme_rid_lookup[gid] if
                        haskey(fluxes_isozymes_backward, rid);
                        init = zero(typeof(enz.value)),
                    ), # flux through negative isozymes
            bound = 0.0,
        ) for (gid, enz) in enzymes if gid in keys(enzyme_rid_lookup)
    )
end

"""
$(TYPEDSIGNATURES)

Helper function to balance the forward or backward isozyme fluxes for a specific
gene product.
"""
function enzyme_balance(
    gid::Symbol,
    rid::Symbol,
    fluxes_isozymes::C.ConstraintTree, # direction
    reaction_isozymes::Dict{String,Dict{String,Isozyme}},
    direction = :kcat_forward,
)
    isozyme_dict = Dict(Symbol(k) => v for (k, v) in reaction_isozymes[string(rid)])

    sum( # this is where the stoichiometry comes in
        -isozyme_value.value *
        isozyme_dict[isozyme_id].gene_product_stoichiometry[string(gid)] /
        getfield(isozyme_dict[isozyme_id], direction) for
        (isozyme_id, isozyme_value) in fluxes_isozymes[rid] if
        gid in Symbol.(keys(isozyme_dict[isozyme_id].gene_product_stoichiometry));
        init = zero(C.LinearValue),
    )
end

function enzyme_capacity(
    enzymes::C.ConstraintTree,
    enzyme_molar_mass::Dict{String,Float64},
    enzyme_ids::Vector{String},
    capacity::Float64,
)
    C.Constraint(
        value = sum(
            enzymes[Symbol(gid)].value * enzyme_molar_mass[gid] for gid in enzyme_ids
        ),
        bound = (0.0, capacity),
    )
end

function build_enzyme_constrained_model(
    model::A.AbstractFBCModel,
    reaction_isozymes::Dict{String,Dict{String,Isozyme}},
    gene_molar_masses::Dict{String,Float64},
    capacity_limitations::Vector{Tuple{String,Vector{String},Float64}},
)
    # create base constraint tree
    m = fbc_model_constraints(model)

    # create directional fluxes
    m +=
        :fluxes_forward^fluxes_in_direction(m.fluxes, :forward) +
        :fluxes_backward^fluxes_in_direction(m.fluxes, :backward)

    # link directional fluxes to original fluxes
    m *=
        :link_flux_directions^sign_split_constraints(
            positive = m.fluxes_forward,
            negative = m.fluxes_backward,
            signed = m.fluxes,
        )

    # create fluxes for each isozyme
    for (rid, _) in m.fluxes_forward
        if haskey(reaction_isozymes, string(rid))
            m +=
                :fluxes_isozymes_forward^rid^isozyme_variables(
                    string(rid),
                    reaction_isozymes,
                )
        end
    end
    for (rid, _) in m.fluxes_backward
        if haskey(reaction_isozymes, string(rid))
            m +=
                :fluxes_isozymes_backward^rid^isozyme_variables(
                    string(rid),
                    reaction_isozymes,
                )
        end
    end

    # link isozyme fluxes to directional fluxes
    m *=
        :link_isozyme_fluxes_forward^link_isozymes(
            m.fluxes_forward,
            m.fluxes_isozymes_forward,
        )
    m *=
        :link_isozyme_fluxes_backward^link_isozymes(
            m.fluxes_backward,
            m.fluxes_isozymes_backward,
        )

    # create enzyme variables
    m += :enzymes^enzyme_variables(model)

    # add enzyme mass balances
    m *=
        :enzyme_stoichiometry^enzyme_stoichiometry(
            m.enzymes,
            m.fluxes_isozymes_forward,
            m.fluxes_isozymes_backward,
            reaction_isozymes,
        )

    # add capacity limitations
    for (id, gids, cap) in capacity_limitations
        m *= Symbol(id)^enzyme_capacity(m.enzymes, gene_molar_masses, gids, cap)
    end

    return m
end
