"""
    FluxVariabilitySummary

A struct used to store summary information about the 
solution of a flux variability analysis result.
"""
struct FluxVariabilitySummary
    biomass_fluxes::Dict{String,Vector{Union{Float64,Nothing}}}
    exchange_fluxes::Dict{String,Vector{Union{Float64,Nothing}}}
end

"""
    flux_variability_summary(flux_result::Tuple{Dict{String, Dict{String, Float64}}, Dict{String, Dict{String, Float64}}}; 
        exclude_exchanges = false,
        exchange_prefixes = _constants.exchange_prefixes,
        biomass_strings = _constants.biomass_strings,
        exclude_biomass = false,
        )::FluxVariabilitySummary

Summarize a dictionary of flux dictionaries obtained eg. from
flux_variability_analysis_dict. The simplified summary representation is useful
for pretty-printing and easily showing the most important results. Internally
this function uses [`looks_like_biomass_reaction`](@ref) and
[`looks_like_exchange_reaction`](@ref). The corresponding keyword arguments
passed to these functions. Use this if your model has non-standard ids for
reactions. 

# Example
```
julia> sol = flux_variability_analysis_dict(model, Gurobi.Optimizer; bounds = objective_bounds(0.99))
julia> flux_res = flux_variability_summary(sol)
Biomass                     Lower bound   Upper bound
  BIOMASS_Ecoli_core_w_GAM: 0.8652        0.8652
Exchange
  EX_h2o_e:                 28.34         28.34
  EX_co2_e:                 22.0377       22.0377
  EX_o2_e:                  -22.1815      -22.1815
  EX_h_e:                   17.3556       17.3556
  EX_glc__D_e:              -10.0         -10.0
  EX_nh4_e:                 -4.8448       -4.8448
  EX_pi_e:                  -3.2149       -3.2149
  EX_for_e:                 0.0           0.0
  ...                       ...           ...
```
"""
function flux_variability_summary(
    flux_result::Tuple{Dict{String,Dict{String,Float64}},Dict{String,Dict{String,Float64}}};
    exclude_exchanges = false,
    exchange_prefixes = _constants.exchange_prefixes,
    biomass_strings = _constants.biomass_strings,
    exclude_biomass = false,
)

    rxn_ids = keys(flux_result[1])
    ex_rxns = filter(
        x -> looks_like_exchange_reaction(
            x,
            exclude_biomass = exclude_biomass,
            biomass_strings = biomass_strings,
            exchange_prefixes = exchange_prefixes,
        ),
        rxn_ids,
    )
    bmasses = filter(
        x -> looks_like_biomass_reaction(
            x;
            exclude_exchanges = exclude_exchanges,
            exchange_prefixes = exchange_prefixes,
            biomass_strings = biomass_strings,
        ),
        rxn_ids,
    )

    biomass_fluxes = Dict{String,Vector{Union{Float64,Nothing}}}()
    for rxn_id in bmasses
        lb = isnothing(flux_result[1][rxn_id]) ? nothing : flux_result[1][rxn_id][rxn_id]
        ub = isnothing(flux_result[2][rxn_id]) ? nothing : flux_result[2][rxn_id][rxn_id]
        biomass_fluxes[rxn_id] = [lb, ub]
    end

    ex_rxn_fluxes = Dict{String,Vector{Union{Float64,Nothing}}}()
    for rxn_id in ex_rxns
        lb = isnothing(flux_result[1][rxn_id]) ? nothing : flux_result[1][rxn_id][rxn_id]
        ub = isnothing(flux_result[2][rxn_id]) ? nothing : flux_result[2][rxn_id][rxn_id]
        ex_rxn_fluxes[rxn_id] = [lb, ub]
    end

    return FluxVariabilitySummary(biomass_fluxes, ex_rxn_fluxes)
end
