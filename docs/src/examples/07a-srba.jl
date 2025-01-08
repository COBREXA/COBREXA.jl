
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
# RBA models require a lot of data. Below we have processed the data into a
# format that is amenable to model building. Sources of the data are shown where
# it is instantiated. In essence, we merely load the information into
# COBREXA/Julia objects.

#md # ```@raw html
#md # <details><summary><strong>Data for reaction turnover numbers</strong></summary>
#md # ```

import CSV
import Statistics: mean

data_root = joinpath(@__DIR__, "src", "examples", "data")

# ### Helper functions
#

# Checks if reaction has a gene reaction rule assigned to it.
function has_reaction_grr(model, rid)
    grr = A.reaction_gene_association_dnf(model, rid)
    if isnothing(grr) || isempty(grr) || isnothing(first(grr)) || isempty(first(grr))
        return false
    end
    return true
end

# This adds kcats to
# 1] transporters (if the kcat is missing and the transporter has a grr)
# 2] physiological reactions which don't have a kcat (but they get an average
#    one if grr is present)
function add_kcats!(
    model,
    kcat_data;
    transporter_kcat = 180.0,
    average_kcat = 25.0,
    block_brackets = false,
)
    for rid in A.reactions(model)
        !has_reaction_grr(model, rid) && continue
        rxn = A.reaction_stoichiometry(model, rid)
        if block_brackets
            suffixes = [split(met, "[")[end] for met in keys(rxn)]
        else
            suffixes = [split(met, "_")[end] for met in keys(rxn)]
        end
        if length(unique(suffixes)) > 1
            if !haskey(kcat_data, rid)
                kcat_data[rid] = transporter_kcat
            end
        end
    end
    for rid in A.reactions(model)
        !has_reaction_grr(model, rid) && continue
        if !haskey(kcat_data, rid)
            kcat_data[rid] = average_kcat
        elseif haskey(kcat_data, rid)
            continue
        end
    end
    return nothing
end

# Finds mass of each gene product in a model
function get_protein_masses(model, proteome_data)
    avg_mass = (mean([x[1] for x in values(proteome_data)]),)
    return Dict{String,Float64}(
        gid => get(proteome_data, gid, avg_mass)[1] for gid in A.genes(model)
    )
end

# Assemble reaction isozymes from kcat data, Uniprot proteome data, and
# ComplexPortal data. This changes the GRR structure of the model (basically
# making it match the expectations from the data files). In a set of complexes,
# this also gets rid of low confidence complexes if any complex in the set can
# be found in the ComplexPortal.
function get_reaction_isozymes!(model, kcat_data, proteome_data, complex_data, scale)
    # protein stoich map, infer from uniprot
    mer_map = Dict(
        "Homomonomer" => 1,
        "Monomer" => 1,
        "Homodimer" => 2,
        "Homotrimer" => 3,
        "Homotetramer" => 4,
        "Homopentamer" => 5,
        "Homohexamer" => 6,
        "Homoheptamer" => 7,
        "Homooctamer" => 8,
        "Homodecamer" => 10,
        "Homododecamer" => 12,
    )
    #+
    # infer protein stoichiometry from uniprot annotations
    reaction_isozymes = Dict{String,Dict{String,Isozyme}}()
    multi_component_enzymes = [] # for use later
    for rid in A.reactions(model)
        haskey(kcat_data, rid) || continue # skip if no kcat data available
        if has_reaction_grr(model, rid)
            grrs = A.reaction_gene_association_dnf(model, rid)
            for (i, grr) in enumerate(grrs)
                isozyme_id = "isozyme_" * string(i)
                d = get!(reaction_isozymes, rid, Dict{String,Isozyme}())
                if length(grr) == 1 # only assign homomers
                    gid = first(grr)
                    if haskey(proteome_data, gid) # has uniprot data
                        idx = proteome_data[gid][2]
                        ss_counts = [get(mer_map, idx, 1.0)]
                    else # no data
                        ss_counts = [1.0]
                    end
                else # assume complexes have uni-stoichiometry, fix in next step
                    ss_counts = fill(1.0, length(grr))
                    push!(multi_component_enzymes, (rid, isozyme_id))
                end
                d[isozyme_id] = Isozyme(
                    gene_product_stoichiometry = Dict(grr .=> ss_counts),
                    kcat_forward = kcat_data[rid] * scale,
                    kcat_reverse = kcat_data[rid] * scale, # assume forward and reverse are the same
                )
            end
        end
    end
    #+
    # fix complex stoichiometry using ComplexPortal data
    fixed_multi_component_enzymes = []
    not_fixed_multi_component_enzymes = []
    for (rid, isozyme_id) in multi_component_enzymes
        grr = collect(keys(reaction_isozymes[rid][isozyme_id].gene_product_stoichiometry))
        stoichs = []
        for v in values(complex_data)
            if length(intersect(collect(keys(v)), grr)) == length(grr) &&
               length(intersect(collect(keys(v)), grr)) == length(v)
                push!(stoichs, v)
            end
        end
        if length(stoichs) == 1
            push!(fixed_multi_component_enzymes, (rid, isozyme_id))
            d = first(stoichs)
            for (k, v) in d
                if v == 0
                    d[k] = 1.0
                end
            end
            reaction_isozymes[rid][isozyme_id].gene_product_stoichiometry = d
        else
            push!(not_fixed_multi_component_enzymes, (rid, isozyme_id))
        end
    end
    #+
    # use only the entries where the stoichiometry is known, if any are known
    for (rid, isozyme_id) in not_fixed_multi_component_enzymes
        if rid in first.(fixed_multi_component_enzymes)
            delete!(reaction_isozymes[rid], isozyme_id)
        end
    end
    return reaction_isozymes
end

# ### Data loading
# load data from papers, complex DB, and Uniprot

proteome_data = begin
    x = CSV.File(joinpath(data_root, "e_coli_uniprot.tsv"), delim = '\t')
    Dict{String,Tuple{Float64,String}}(
        gp => (m, coalesce(k, "")) for (gp, m, k) in zip(x.gene_product, x.mass, x.kind)
    )
end

# this data is taken from: *Heckmann, David, et al. "Machine learning applied
# to enzyme turnover numbers reveals protein structural correlates and improves
# metabolic models." Nature communications 9.1 (2018): 1-10.*
kcat_data = Dict(
    String(r.KcatID) => r.KcatHeckmann for
    r in CSV.File(joinpath(data_root, "ecoli_kcats.tsv"))
)

complex_data = begin
    x = CSV.File(joinpath(data_root, "e_coli_complex.tsv"), delim = '\t')
    res = Dict{String,Dict{String,Float64}}()
    for c in x.complex
        res[c] = Dict{String,Float64}()
    end
    for (c, gp, s) in zip(x.complex, x.gene_product, x.stoichiometry)
        res[c][gp] = s
    end
    res
end

# load the E. coli model as a basis to assemble the data around

model = load_model("iML1515.json", A.CanonicalModel.Model)

# Assemble the isozymes and gene product masses increase protein concentrations
# by changing units, convert kcat units from 1/s to k/h
scale = 3600 / 1e3

# add kcats that are not measured
# https://bionumbers.hms.harvard.edu/bionumber.aspx?id=114686&ver=3&trm=glucose+transport&org=
add_kcats!(model, kcat_data; transporter_kcat = 180.0, average_kcat = 25.0)

# get uniprot molar masses (units: kDa = kg/mol)
gene_product_molar_masses = get_protein_masses(model, proteome_data)

# collect kcats and stoichiometry
reaction_isozymes =
    get_reaction_isozymes!(model, kcat_data, proteome_data, complex_data, scale)

# apply analysis-specific quirks
gene_product_molar_masses["b1692"] = 54.58
gene_product_molar_masses["s0001"] = 54.58

# ribosome data (see comment URLs for sources)
# https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100118&ver=10&trm=e+coli+ribosome+molar+mass&org= )
gene_product_molar_masses["ribosome"] = 2700.0
# https://bionumbers.hms.harvard.edu/bionumber.aspx?id=101175
aas_in_ribosome = 7459
# https://bionumbers.hms.harvard.edu/bionumber.aspx?id=114971&ver=1&trm=ribosome+amino+acid&org=
atp_polymerization_cost = 4.2
# https://bionumbers.hms.harvard.edu/bionumber.aspx?id=109522&ver=2&trm=ribosome+density&org=
ribosomes_per_protein = 3.46 # multiple ribosomes can produce a single AA, effectively speeding up the overall rate (assume average for all genes)

# Load and prepare the amino acid data

amino_acids = Dict(
    :glu__L_c => "E",
    :lys__L_c => "K",
    :his__L_c => "H",
    :ser__L_c => "S",
    :arg__L_c => "R",
    :gln__L_c => "Q",
    :trp__L_c => "W",
    :ile__L_c => "I",
    :tyr__L_c => "Y",
    :gly_c => "G",
    :thr__L_c => "T",
    :ala__L_c => "A",
    :cys__L_c => "C",
    :asn__L_c => "N",
    :pro__L_c => "P",
    :phe__L_c => "F",
    :asp__L_c => "D",
    :val__L_c => "V",
    :leu__L_c => "L",
    :met__L_c => "M",
)

aacount = begin # count of each amino acid in a gene product
    x = CSV.File(joinpath(data_root, "ecoli_gene_product_aa_counts.tsv"), delim = '\t')
    res = Dict{Symbol,Dict{Symbol,Int}}()
    for gp in x.gene_product
        res[Symbol(gp)] = Dict{Symbol,Int}()
    end
    for (gp, met, n) in zip(x.gene_product, x.metabolite, x.aa_count)
        res[Symbol(gp)][Symbol(met)] = n
    end
    res
end

avg_prot = Dict()
for v in values(aacount)
    for (k, vv) in v
        avg_prot[k] = get(avg_prot, k, 0) + vv
    end
end
for (k, v) in avg_prot
    avg_prot[k] = v / sum(values(avg_prot))
end
aacount[:ribosome] =
    Dict(k => round(Int, aas_in_ribosome * v, RoundNearest) for (k, v) in avg_prot)

#md # ```@raw html
#md # </details>
#md # ```

# ## Model assembly
# Now that we have loaded all the data, we can start building the simplified RBA model

# First, load the model
model = load_model("iML1515.json", A.CanonicalModel.Model)

# Next, we will remove the biomass reactions, as these are handled differently
# in RBA. However, we will save biomass composition for later use.
biomass = model.reactions["BIOMASS_Ec_iML1515_core_75p37M"]
# remove the biomass functions for rba
delete!(model.reactions, "BIOMASS_Ec_iML1515_WT_75p37M")
delete!(model.reactions, "BIOMASS_Ec_iML1515_core_75p37M")

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
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = [
        ("membrane", membrane_gids, float(total_capacity * membrane_frac)),
        ("total", A.genes(model), float(total_capacity)),
    ],
)

# apply a few quirks
ct.fluxes[:EX_glc__D_e].bound = C.Between(-1000, 1000)
ct.fluxes_forward[:EX_glc__D_e].bound = C.Between(0, 1000)
ct.fluxes_reverse[:EX_glc__D_e].bound = C.Between(0, 1000)
ct.fluxes[:EX_for_e].bound = C.EqualTo(0.0)
ct.fluxes[:EX_pyr_e].bound = C.EqualTo(0.0)
ct.fluxes[:EX_5dglcn_e].bound = C.EqualTo(0.0)
ct.fluxes[:EX_lac__D_e].bound = C.EqualTo(0.0)

# We now need to extend the enzyme constrained model with further resource
# allocation constraints. Since RBA models are bilinear, it is simpler to create
# a function that builds a model at a specific growth rate, transforming the
# problem into a linear one. Future extensions of cobrexa will allow parameters
# to be inserted into models directly, obviating the need for this
# function-based appraoch. The function belows glues the sRBA constraint system
# to the given ecRBA model at a specific growth rate.
function with_srba_constraints(ct, mu)
    #+
    # atp + h2o -> adp + h + po4
    energy_metabolites =
        Dict(:atp_c => -1, :h2o_c => -1, :adp_c => 1, :h_c => 1, :pi_c => 1)
    #+
    # attach new variables for ribosomes (recall ribosomes are required to make proteins, and ribosomes themselves)
    rbatree =
        ct +
        :ribosomes^C.variables(;
            keys = [collect(keys(ct.gene_product_amounts)); :ribosome],
            bounds = C.Between(0, Inf),
        )
    #+
    # we are going to modify a few parts locally so let's make a copy of them to not change the base model (ct)
    rbatree = C.ConstraintTree(
        rbatree...,
        :flux_stoichiometry => deepcopy(rbatree.flux_stoichiometry),
        :gene_product_capacity => deepcopy(rbatree.gene_product_capacity),
    )
    #+
    # amino acid dilution by growth note: in RBA biomass components get diluted
    # directly, and it is not necessary to have a biomass reaction for this
    # purpose. Here the amino acids get special attention, since they are
    # consumed by the ribosomes.
    for aa in keys(amino_acids)
        rbatree.flux_stoichiometry[aa].value -=
            mu *
            0.001 *
            (
                sum(
                    aacount[g][aa] * c.value for (g, c) in rbatree.gene_product_amounts if
                    haskey(aacount, g) && haskey(aacount[g], aa);
                    init = zero(C.LinearValue),
                ) + sum(
                    aacount["ribosome"][String(aa)] * C.value(c) for
                    (_, c) in rbatree.ribosomes if haskey(aacount[:ribosome], String(aa));
                    init = zero(C.LinearValue),
                )
            )
    end
    #+
    # In contrast, the other biomass components just get diluted by growth, including the energy metabolites
    prot_atp_req = 12.0 # protein polymerization cost for average proteome in 1 gDW
    for (mid, v) in biomass.stoichiometry # from the saved biomass composition of the original model
        k = Symbol(mid)
        haskey(amino_acids, k) && continue
        #+
        # if the metabolite is an energy metabolite, we add the energy cost of translation
        offset = haskey(energy_metabolites, k) ? -energy_metabolites[k] * prot_atp_req : 0.0
        rbatree.flux_stoichiometry[k].value += mu * (v + offset)
    end
    #+
    # energy costs of protein synthesis
    for (k, kk) in energy_metabolites
        rbatree.flux_stoichiometry[k].value +=
            (kk * mu * 0.001 * atp_polymerization_cost) * sum(
                sum(values(aacount[g])) * C.value(c) for
                (g, c) in rbatree.gene_product_amounts if haskey(aacount, g);
                init = C.zero(C.LinearValue),
            )
        #+
        # ribosome synthesis also costs energy
        rbatree.flux_stoichiometry[k].value +=
            (kk * mu * 0.001 * atp_polymerization_cost * aas_in_ribosome) *
            sum(C.value(x) for x in values(rbatree.ribosomes); init = C.zero(C.LinearValue))
    end
    #+
    # add the equations for protein synthesis by ribosomes
    # (the conversion is: ribosome elongation rate amino acids/second => amino acids/hr)
    kr = ribosomes_per_protein * 12 * 3600
    aa_sum(g) = haskey(aacount, g) ? sum(values(aacount[g])) : 300
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
    # patch the total gene product capacity bound with the new ribosomes
    rbatree.gene_product_capacity.total.value +=
        sum(C.value.(values(rbatree.ribosomes))) * gene_product_molar_masses["ribosome"]

    return rbatree
end


# ## Run the simulations
# Here we use screen to efficiently run all the simulations through the
# expectably viable growth range

mus = range(0.1, 0.95, 10) # simulate at these growth rates

settings = [
    set_optimizer_attribute("solver", "simplex"),
    set_optimizer_attribute("simplex_strategy", 4),
]

res = screen(mus) do mu
    @info "sRBA step" mu
    rbat = with_srba_constraints(ct, mu)
    rbat *=
        :lp_objective^C.Constraint(
            sum(
                C.value(v) * gene_product_molar_masses[string(k)] for
                (k, v) in rbat.gene_product_amounts
            ) + sum(
                C.value(v) * gene_product_molar_masses["ribosome"] for
                v in values(rbat.ribosomes)
            ),
            nothing,
        )
    sol = optimized_values(
        rbat;
        settings = [silence],
        objective = rbat.lp_objective.value,
        optimizer = HiGHS.Optimizer,
        sense = Minimal,
        settings,
    )
    isnothing(sol) && return nothing
    return (;
        mu,
        ribosome_mass = gene_product_molar_masses["ribosome"] * sum(values(sol.ribosomes)),
        total_mass = sol.gene_product_capacity.total,
        ac_flux = sol.fluxes.EX_ac_e,
        glc_flux = sol.fluxes.EX_glc__D_e,
        o2_flux = sol.fluxes.EX_o2_e,
    )
end

# finally, we can plot the data, to see if we can recapitulate known phenomena

#= TODO

using CairoMakie

# load measured ribosome protein mass fractions
ribosome_measurements = CSV.File(joinpath(data_root, "ecoli_ribosomes.tsv"))

# First, show that the predicted ribosome density matches experimental
# observations, and also show that overflow metabolism occurs (latter is due to
# the membrane bound).
fig = Figure();
ax = Axis(fig[1, 1], xlabel = "Growth rate, 1/h", ylabel = "Ribosome mass fraction")
scatter!(
    ax,
    ribosome_measurements.GR,
    ribosome_measurements.ActiveR ./ 100,
    label = "Measurements",
)
lines!(ax, mus, [r.ribosome_mass / r.total_mass for r in res], label = "Model predictions")
axislegend(ax, position = :lt)

ax2 = Axis(fig[2, 1], xlabel = "Growth rate, 1/h", ylabel = "Metabolite flux")
lines!(ax2, mus, [r.ac_flux for r in res], label = "Acetate")
lines!(ax2, mus, [abs(r.glc_flux) for r in res], label = "Glucose")
lines!(ax2, mus, [abs(r.o2_flux) for r in res], label = "Oxygen")
axislegend(ax2, position = :lt)
fig

=#
