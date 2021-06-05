"""
    atom_exchange(flux_dict::Dict{String, Float64}, model::StandardModel)

Return a dictionary mapping the flux of atoms across the boundary of the model 
given `flux_dict` (the solution of a constraint based analysis) of reactions in `model`.
"""
function atom_exchange(flux_dict::Dict{String,Float64}, model::StandardModel)
    atom_flux = Dict{String,Float64}()
    for (rxn_id, flux) in flux_dict
        if is_boundary(model.reactions[rxn_id])
            for (met, stoich) in model.reactions[rxn_id].metabolites
                adict = get_atoms(model.metabolites[met])
                for (atom, stoich) in adict
                    atom_flux[atom] = get(atom_flux, atom, 0.0) + flux * stoich
                end
            end
        end
    end
    return atom_flux
end

"""
    metabolite_fluxes(flux_dict::Dict{String, Float64}, model::StandardModel)

Return two dictionaries of metabolite `id`s mapped to reactions that consume or 
produce them given the flux distribution supplied in `fluxdict`.
"""
function metabolite_fluxes(flux_dict::Dict{String,Float64}, model::StandardModel)
    S = stoichiometry(model)
    met_flux = Dict{String,Float64}()
    rxnids = reactions(model)
    metids = metabolites(model)

    producing = Dict{String,Dict{String,Float64}}()
    consuming = Dict{String,Dict{String,Float64}}()
    for (row, metid) in enumerate(metids)
        for (col, rxnid) in enumerate(rxnids)
            mf = flux_dict[rxnid] * S[row, col]
            # ignore zero flux
            if mf < -_constants.tolerance # consuming rxn
                if haskey(consuming, metid)
                    consuming[metid][rxnid] = mf
                else
                    consuming[metid] = Dict(rxnid => mf)
                end
            elseif mf > _constants.tolerance
                if haskey(producing, metid)
                    producing[metid][rxnid] = mf
                else
                    producing[metid] = Dict(rxnid => mf)
                end
            end
        end
    end
    return consuming, producing
end

"""
    set_bound(index, optimization_model;
        ub=_constants.default_reaction_rate,
        lb=-_constants.default_reaction_rate)

Helper function to set the bounds of variables.
The JuMP `set_normalized_rhs` function is a little confusing, 
so this function simplifies setting constraints. In short, JuMP
uses a normalized right hand side representation of constraints, 
which means that lower bounds have their sign flipped. This function
does this for you, so you don't have to remember to do this whenever you
change the constraints. 

Just supply the constraint `index` and the JuMP model (`opt_model`) that 
will be solved, and the variable's bounds will be set to `ub` and `lb`.
"""
function set_bound(
    vind,
    opt_model;
    ub = _constants.default_reaction_rate,
    lb = -_constants.default_reaction_rate,
)
    set_normalized_rhs(opt_model[:lbs][vind], -lb)
    set_normalized_rhs(opt_model[:ubs][vind], ub)
end
