
using COBREXA, JSON, GLPK
using JSONFBCModels
import AbstractFBCModels as A
import ConstraintTrees as C

using Gurobi

# get model
download_model(
    "http://bigg.ucsd.edu/static/models/iML1515.json",
    "iML1515.json",
    "b0f9199f048779bb08a14dfa6c09ec56d35b8750d2f99681980d0f098355fbf5",
)

model = load_model("iML1515.json")

rid_kcats = JSON.parsefile(joinpath("docs", "src", "showcase", "data", "ecoli_kcats.json"))
reaction_isozymes = Dict{String,Dict{String,Isozyme}}() # a mapping from reaction IDs to isozyme IDs to isozyme structs.
for rid in A.reactions(model)
    grrs = A.reaction_gene_association_dnf(model, rid)
    isnothing(grrs) && continue # skip if no grr available
    for (i, grr) in enumerate(grrs)
        d = get!(reaction_isozymes, rid, Dict{String,Isozyme}())
        kcat_f = get(rid_kcats, rid, 30.0)
        kcat_r = get(rid_kcats, rid * "_b", kcat_f)
        d["isozyme_"*string(i)] = Isozyme(
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), # assume subunit stoichiometry of 1 for all isozymes
            kcat_forward = kcat_f * 3.6, # forward reaction turnover number units = k/h
            kcat_reverse = kcat_r * 3.6, # reverse reaction turnover number units = k/h
        )
    end
end

gene_product_molar_masses = Dict(
    k => v for
    (k, v) in JSON.parsefile(joinpath("docs", "src", "showcase", "data", "ecoli_gpmms.json"))
)
gene_product_molar_masses["b1692"] = 54.58 # missing
gene_product_molar_masses["s0001"] = 54.58 # missing

ecoli_measured_abundances = JSON.parsefile(joinpath("docs", "src", "showcase", "data", "ecoli_community_abundances.json"))

# using smoment, so just KO the reaction
kos = Dict(
    "argA" => (:ACGS, :arg),
    "serA" => (:PGCD, :ser),
    "ilvA" => (:THRD_L, :ile),
    "metA" => (:HSST, :met),
    "hisB" => (:HISTP, :his),
    "proA" => (:G5SD, :pro),
    "cysE" => (:SERAT, :cys),
    "thrC" => (:THRS, :thr),
    "glyA" => (:GHMT2r, :gly),
    "leuB" => (:IPMD, :leu),
    "trpC" => (:IGPS, :trp),
    "pheA" => (:PPNDH, :phe),
    "tyrA" => (:PPND, :tyr),
    "lysA" => (:DAPDC, :lys),
)

ko_rxn1, aa1 = kos["ilvA"]
ko_rxn2, aa2 = kos["metA"]

aa1_id = aa1 == :gly ? Symbol("EX_", aa1, "_e") : Symbol("EX_", aa1, "__L_e")
aa2_id = aa2 == :gly ? Symbol("EX_", aa2, "_e") : Symbol("EX_", aa2, "__L_e")

open_bounds!(model, rxn_id) = begin
    model.fluxes[rxn_id].bound = C.Between(-1000,1000)
    model.fluxes_forward[rxn_id].bound = C.Between(0,1000) # TODO unfortunate, easy to miss
    model.fluxes_reverse[rxn_id].bound = C.Between(0,1000) # TODO unfortunate, easy to miss
end

mut1 = simplified_enzyme_constrained_flux_balance_constraints(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = 400.0, # mg/gDW
    interface = :identifier_prefixes,
)
open_bounds!(mut1, :EX_glc__D_e)
open_bounds!(mut1, aa1_id)
open_bounds!(mut1, aa2_id)

mut2 = simplified_enzyme_constrained_flux_balance_constraints(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = 400.0, # mg/gDW
    interface = :identifier_prefixes,
)
open_bounds!(mut2, :EX_glc__D_e)
open_bounds!(mut2, aa1_id)
open_bounds!(mut2, aa2_id)

# create bound function to constrain interface
boundf(id) = begin
    ex_id = first(id)
    if ex_id == aa1_id || ex_id == aa2_id
        C.EqualTo(0.0)
    else
        mut1.interface.exchanges[ex_id].bound # have same interface, so easy
    end
end

a = 0.80
x = interface_constraints(
    :mut1 => (mut1, mut1.interface.exchanges, a),
    :mut2 => (mut2, mut2.interface.exchanges, 1-a);
    bound = boundf
)

x *= :equalgrowth^equal_value_constraint(x.mut1.fluxes.BIOMASS_Ec_iML1515_core_75p37M, x.mut2.fluxes.BIOMASS_Ec_iML1515_core_75p37M)
    
sol = optimized_values(
    x;
    optimizer = Gurobi.Optimizer,
    objective = x.mut1.objective.value,
    sense = Maximal,
    settings = [silence,]
)
sol.mut1.objective
