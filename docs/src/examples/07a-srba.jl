
# Copyright (c) 2021-2024, University of Luxembourg                         #src
# Copyright (c) 2021-2024, Heinrich-Heine University Duesseldorf            #src
#                                                                           #src
# Licensed under the Apache License, Version 2.0 (the "License");           #src
# you may not use this file except in compliance with the License.          #src
# You may obtain a copy of the License at                                   #src
#                                                                           #src
#     http://www.apache.org/licenses/LICENSE-2.0                            #src
#                                                                           #src
# Unless required by applicable law or agreed to in writing, software       #src
# distributed under the License is distributed on an "AS IS" BASIS,         #src
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  #src
# See the License for the specific language governing permissions and       #src
# limitations under the License.                                            #src

# # Simplified resource balance analysis
# Resource balance analysis models are an extension of enzyme constrained
# models. They can incorporate translation, transcription and replication
# effects to dramatically increase the realism of a metabolic model. To
# demonstrate the ease of incrementally building complex models in COBREXA, we
# will build a simplified resource balance analysis model from the ground up,
# primarily focussing on translation extensions to an enzyme constrained flux
# balance analysis model.

using COBREXA

# In contrast to the previous examples, we will use the iML1515 model of _E.
# coli_, which is the latest, full-scale genome-scale metabolic model of the
# organism.

download_model(
    "http://bigg.ucsd.edu/static/models/iML1515.json",
    "iML1515.json",
    "b0f9199f048779bb08a14dfa6c09ec56d35b8750d2f99681980d0f098355fbf5",
)

# In additional to COBREXA, the optimization solver, and the model format
# package, we will use ConstraintTrees to simplify the RBA model construction
# process:

import AbstractFBCModels as A
import JSONFBCModels
import ConstraintTrees as C
import HiGHS

# ## Collect data for RBA model

import CSV
import Statistics: mean

data_dir = joinpath(dirname(pathof(COBREXA)), "..", "docs", "src", "examples", "data")

e_coli_gp_mass = Dict{String,Float64}(
    x.gene_product => x.mass for
    x in CSV.File(joinpath(data_dir, "e_coli_gp_mass.tsv"), delim = '\t')
)

kcat_scale = 3600 / 1e3
e_coli_rxn_kcat_isozyme = Dict{String,Isozyme}(
    x.reaction => Isozyme(
        kcat_forward = x.kcat * kcat_scale,
        kcat_reverse = x.kcat * kcat_scale,
        gene_product_stoichiometry = Dict(),
    ) for x in CSV.File(joinpath(data_dir, "e_coli_reaction_kcat.tsv"), delim = '\t')
)

e_coli_rxn_isozymes = Dict{String,Dict{String,Isozyme}}()
for x in CSV.File(joinpath(data_dir, "e_coli_isozyme_gp.tsv"), delim = '\t')
    haskey(e_coli_rxn_kcat_isozyme, x.reaction) || continue
    rxn = get!(e_coli_rxn_isozymes, x.reaction, Dict{String,Isozyme}())
    iso = get!(rxn, x.isozyme, deepcopy(e_coli_rxn_kcat_isozyme[x.reaction]))
    iso.gene_product_stoichiometry[x.gene_product] = x.stoichiometry
end

e_coli_gp_aas = Dict{String,Dict{Symbol,Int}}(
    begin
        d = Dict(keys(x) .=> values(x))
        gp = d[:gene_product]
        delete!(d, :gene_product)
        gp => d
    end for x in CSV.File(joinpath(data_dir, "e_coli_gp_aa.tsv"), delim = '\t')
)

amino_acids = Set(aa for (k, v) in e_coli_gp_aas for (aa, _) in v)

# assume some constants (these can be found via https://bionumbers.hms.harvard.edu )
aas_in_ribosome = sum(values(e_coli_gp_aas["ribosome"]))
atp_polymerization_cost = 4.2
protein_polymerization_atp_per_gDW = 12.0

# ## Model assembly

# First, load the model
model = load_model("iML1515.json", A.CanonicalModel.Model)

# We remove the biomass reactions, as these are handled differently in RBA.
# However, we will save biomass composition for later use.
biomass = model.reactions["BIOMASS_Ec_iML1515_core_75p37M"]
# remove the biomass functions for rba
delete!(model.reactions, "BIOMASS_Ec_iML1515_WT_75p37M")
delete!(model.reactions, "BIOMASS_Ec_iML1515_core_75p37M")

# apply a few quirks to enforce glucose-consuming metabolism
model.reactions["EX_glc__D_e"].lower_bound = -1000
model.reactions["EX_glc__D_e"].upper_bound = 1000
delete!.(Ref(model.reactions), ["EX_for_e", "EX_pyr_e", "EX_5dglcn_e", "EX_lac__D_e"]);

# identify membrane reactions (transporting stuff between compartments)
membrane_rids = [
    rid for (rid, r) in model.reactions if
    length(unique(last.(split.(keys(r.stoichiometry), "_")))) != 1
]

# identify membrane proteins as the ones catalyzing the membrane reactions from above
membrane_gids = unique(
    g for
    rid in membrane_rids if !isnothing(A.reaction_gene_association_dnf(model, rid)) for
    gs in A.reaction_gene_association_dnf(model, rid) for g in gs
)

# Now we build an enzyme constrained metabolic model using COBREXA. Note, we
# assume two capacity limitations here: 1) a total capacity bound, and a
# membrane capacity bound.
total_capacity = 550.0
membrane_frac = 0.20 # fraction of total proteome that is allowed to be membrane associated

# use COBREXA's built in function to lay the foundation of the RBA model
ct = enzyme_constrained_flux_balance_constraints(
    model;
    reaction_isozymes = e_coli_rxn_isozymes,
    gene_product_molar_masses = e_coli_gp_mass,
    capacity = Tuple{String,Vector{String},Float64}[
        ("membrane", membrane_gids, float(total_capacity * membrane_frac)),
        ("total", A.genes(model), float(total_capacity)),
    ],
)

# How ATP gets consumed (this is used to power the RBA processes).
energy_stoichiometry = Dict(:atp_c => -1, :h2o_c => -1, :adp_c => 1, :h_c => 1, :pi_c => 1)

# We now need to extend the enzyme constrained model with further resource
# allocation constraints. Since RBA models are bilinear, it is simpler to create
# a function that builds a model at a specific growth rate, transforming the
# problem into a linear one. Future extensions of cobrexa will allow parameters
# to be inserted into models directly, obviating the need for this
# function-based appraoch. The function belows glues the sRBA constraint system
# to the given ecRBA model at a specific growth rate.
function with_srba_constraints(ct, mu; ribosome_aa_per_second = 12)
    #+
    # First, attach new variables for ribosome production (recall ribosomes are
    # required to make proteins, and ribosomes themselves)
    rbatree =
        ct +
        :ribosomes^C.variables(;
            keys = [collect(keys(ct.gene_product_amounts)); :ribosome],
            bounds = C.Between(0, Inf),
        )
    #+
    # we are going to modify a few parts locally so let's make a copy of them
    # to not change the base model (ct)
    rbatree = C.ConstraintTree(
        rbatree...,
        :flux_stoichiometry => deepcopy(rbatree.flux_stoichiometry),
        :gene_product_capacity => deepcopy(rbatree.gene_product_capacity),
    )
    #+
    # Connect amino-acid production to its consumption for gene products.
    #
    # (The conversion factor 0.001 is required to convert the units from mol to
    # mmol, which is assumed in AA-related constraints.)
    for aa in amino_acids
        rbatree.flux_stoichiometry[aa].value -=
            mu *
            0.001 *
            (
                sum(
                    e_coli_gp_aas[g][aa] * c.value for
                    (g, c) in rbatree.gene_product_amounts if
                    haskey(e_coli_gp_aas, g) && haskey(e_coli_gp_aas[g], aa);
                    init = zero(C.LinearValue),
                ) + sum(
                    e_coli_gp_aas["ribosome"][String(aa)] * C.value(c) for
                    (_, c) in rbatree.ribosomes if
                    haskey(e_coli_gp_aas["ribosome"], String(aa));
                    init = zero(C.LinearValue),
                )
            )
    end
    #+
    # Make a modified biomass reaction which is adjusted for growth (mu),
    # and additionally consumes the energy metabolites required for the
    # polymerization.
    for (mid, v) in biomass.stoichiometry
        k = Symbol(mid)
        k in amino_acids && continue
        offset =
            haskey(energy_stoichiometry, k) ?
            -energy_stoichiometry[k] * protein_polymerization_atp_per_gDW : 0.0
        rbatree.flux_stoichiometry[k].value += mu * (v + offset)
    end
    #+
    # Consume energy (ATP) for all polymerization.
    for (k, kk) in energy_stoichiometry
        rbatree.flux_stoichiometry[k].value +=
            (kk * mu * 0.001 * atp_polymerization_cost) * sum(
                sum(values(e_coli_gp_aas[g])) * C.value(c) for
                (g, c) in rbatree.gene_product_amounts if haskey(e_coli_gp_aas, g);
                init = C.zero(C.LinearValue),
            )
        #+
        # ribosome synthesis also costs energy
        rbatree.flux_stoichiometry[k].value +=
            (kk * mu * 0.001 * atp_polymerization_cost * aas_in_ribosome) *
            sum(C.value(x) for x in values(rbatree.ribosomes); init = C.zero(C.LinearValue))
    end
    #+
    # Use and consume ribosomes for protein synthesis.
    # (the conversion below is: ribosome elongation rate amino acids/second => amino acids/hr)
    kr = ribosome_aa_per_second * 3600
    aa_sum(g) = haskey(e_coli_gp_aas, g) ? sum(values(e_coli_gp_aas[g])) : 300
    rbatree *=
        :protein_synthesis^C.ConstraintTree(
            g =>
                C.Constraint(kr * rbatree.ribosomes[g].value - mu * c.value * aa_sum(g), 0)
            for (g, c) in rbatree.gene_product_amounts
        )
    #+
    # add the equation for ribosome synthesis by themselves
    rbatree *=
        :ribosome_balance^C.Constraint(
            kr * rbatree.ribosomes[:ribosome].value -
            (mu * aas_in_ribosome) *
            sum(C.value.(values(rbatree.ribosomes)); init = zero(C.LinearValue)),
            0,
        )
    #+
    # Compute the total ribosome production mass (which can
    # serve as minimization objective)
    rbatree *=
        :total_ribosome_mass^C.Constraint(
            C.sum(
                (
                    C.value(v) * e_coli_gp_mass["ribosome"] for
                    v in values(rbatree.ribosomes)
                ),
                init = zero(C.LinearValue),
            ),
        )
    #+
    # Add the ribosome mass into the total capacity bound
    rbatree.gene_product_capacity.total.value += rbatree.total_ribosome_mass.value
    return rbatree
end


# ## Run the simulations
# Here we use screen to efficiently run all the simulations through the
# expectably viable growth range

mus = range(0.1, 1.2, 10) # simulate at these growth rates

@time res = screen(mus, workers = [1]) do mu
    rba_constraints = with_srba_constraints(ct, mu, ribosome_aa_per_second = 12)
    sol = optimized_values(
        rba_constraints;
        objective = rba_constraints.gene_product_capacity.total.value,
        sense = Minimal,
        optimizer = HiGHS.Optimizer,
    )
    isnothing(sol) && return nothing
    return (;
        mu,
        membrane_mass = sol.gene_product_capacity.membrane,
        ribosome_mass = sol.total_ribosome_mass,
        enzyme_mass = sol.gene_product_capacity.total - sol.total_ribosome_mass,
        total_mass = sol.gene_product_capacity.total,
        ac_flux = sol.fluxes.EX_ac_e,
        glc_flux = sol.fluxes.EX_glc__D_e,
        o2_flux = sol.fluxes.EX_o2_e,
    )
end
