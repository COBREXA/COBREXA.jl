
using COBREXA, JSON, GLPK
using JSONFBCModels
import AbstractFBCModels as A

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

mut1 = simplified_enzyme_constrained_flux_balance_constraints(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = 400.0, # mg/gDW
)
mut1.fluxes.EX_glc__D_e.bound = C.Between(-1000, 1000)
mut1.fluxes[ko_rxn1].bound = C.EqualTo(0.0) # knockout
mut1.fluxes[aa1_id].bound = C.Between(-1000,1000)
mut1.fluxes[aa2_id].bound = C.Between(-1000,1000)

