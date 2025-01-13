
# Copyright (c) 2021-2025, University of Luxembourg                         #src
# Copyright (c) 2021-2025, Heinrich-Heine University Duesseldorf            #src
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

# # Enzyme constrained community flux balance analysis
# In this example we will demonstrate how to merge two enzyme constrained models
# together, and run community flux balance analysis simulation on them. The goal
# is to compare the optimal community abundance of models using cFBA and ec-cFBA
# simulations, illustrating the superiority of the latter approach. The specific
# communities simulated here will be _E.coli_ auxotrophic co-cultures. Each
# community member is dependent on the other through an amino acid knockout. The
# simulations are based on the experimental work in *Mee, Michael T., et al.
# "Syntrophic exchange in synthetic microbial communities." Proceedings of the
# National Academy of Sciences 111.20 (2014)*. Finally, the predicted abundances
# using the two modeling approaches will be compared to experimental data.

using COBREXA

# Like the previous tutorial, we will use the iML1515 model of _E.
# coli_, which is the latest, full-scale genome-scale metabolic model of the
# organism.

download_model(
    "http://bigg.ucsd.edu/static/models/iML1515.json",
    "iML1515.json",
    "b0f9199f048779bb08a14dfa6c09ec56d35b8750d2f99681980d0f098355fbf5",
)

# In additional to COBREXA, the optimization solver, and the model format
# package, we will use ConstraintTrees to simplify the model construction
# process:

import AbstractFBCModels as A
import JSONFBCModels
import ConstraintTrees as C
import HiGHS

# ## Collect data for enzyme constrained models
# Like in the previous example, we will spend some time gathering the necessary
# data to build the wild type enzyme constrained metabolic models.

#md # ```@raw html
#md # <details><summary><strong>Data for reaction turnover numbers</strong></summary>
#md # ```

import CSV
import Statistics: mean

data_root = joinpath(dirname(pathof(COBREXA)), "..", "docs", "src", "examples", "data")

# ### Helper functions

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
                # println("Transporter ", rid, " assigned kcat.")
                kcat_data[rid] = transporter_kcat
            end
        end
    end

    for rid in A.reactions(model)
        !has_reaction_grr(model, rid) && continue
        if !haskey(kcat_data, rid)
            kcat_data[rid] = average_kcat
            # println("Assigned average kcat to: ", rid)
        elseif haskey(kcat_data, rid)
            continue
        else
            # println("Rid not assigned a kcat: ", rid)
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

#md # ```@raw html
#md # </details>
#md # ```

# ## Build community models
# Now we will focus on building the community models. Since constraints can be
# added incrementally, and are not dependent on the modeling framework, we only
# need to write a general purpose community construction function, which will
# naturally work with either classic FBA or enzyme constrained FBA models.

# below we create a community assembly function, assuming an input of a wildtype model
function auxotrophe_community(wt, aa_ko; fbc_only = false)

    # these are the possible KOs we consider
    kos = Dict(
        :argA => (:ACGS, :arg),
        :ilvA => (:THRD_L, :ile),
        :metA => (:HSST, :met),
        :hisB => (:HISTP, :his),
        :proA => (:G5SD, :pro),
        :cysE => (:SERAT, :cys),
        :thrC => (:THRS, :thr),
        :leuB => (:IPMD, :leu),
        :trpC => (:IGPS, :trp),
        :pheA => (:PPNDH, :phe),
        :tyrA => (:PPND, :tyr),
        :lysA => (:DAPDC, :lys),
    )

    aa_id(aa) = Symbol("EX_", aa, "__L_e")

    # bound the exchanges
    boundf(id) = begin
        ex_id = first(id)
        if fbc_only && ex_id == :EX_glc__D_e
            C.Between(-50, 0)
        else
            wt.interface.exchanges[ex_id].bound
        end
    end

    # create KO pairs to be joined via their exchanges
    ko_pairs = [
        ko => (deepcopy(wt), deepcopy(wt.interface.exchanges), ko_ab) for
        (ko, ko_ab) in aa_ko
    ]

    # create community model (interface joins them)
    x = interface_constraints(ko_pairs...; bound = boundf)

    all_kos = first.(aa_ko) # get model ids
    x *= # set all growth rates equal
        :equalgrowth^all_equal_constraints(
            x[first(all_kos)].fluxes.BIOMASS_Ec_iML1515_core_75p37M,
            C.ConstraintTree(
                ko => x[ko].fluxes.BIOMASS_Ec_iML1515_core_75p37M for ko in all_kos[2:end]
            ),
        )
    x *= # set objective
        :objective^C.Constraint(
            x[first(all_kos)].fluxes.BIOMASS_Ec_iML1515_core_75p37M.value,
            nothing,
        )

    # knockout - single gene
    for ko in all_kos
        x[ko].fluxes[first(kos[ko])].bound = C.EqualTo(0)
    end

    # allow exchanges of KO'd aa
    all_aa_exs = [aa_id(last(kos[ko])) for ko in all_kos]
    for ko in all_kos
        for aa_ex in all_aa_exs
            if fbc_only
                open_bounds_fbc!(x[ko], aa_ex)
            else
                open_bounds_ecfbc!(x[ko], aa_ex)
            end
        end
    end

    x
end

# simulate a model
function auxotrophe_fba(wt, aa_ko; fbc_only = false)
    x = auxotrophe_community(wt, aa_ko; fbc_only)

    sol = optimized_values(
        x;
        optimizer = HiGHS.Optimizer,
        objective = x.objective.value,
        sense = Maximal,
    )
    isnothing(sol) && return nothing

    (; mu = sol.objective, (Symbol(aa) => ab for (aa, ab) in aa_ko)...)
end

# open specific bounds
function open_bounds_fbc!(ct, rxn_id)
    ct.fluxes[rxn_id].bound = C.Between(-1000, 1000)
    ct.interface.exchanges[rxn_id].bound = C.Between(-1000, 1000)
end

# since enzyme constrained models add more reactions to the models, the bound opening is more intensive
function open_bounds_ecfbc!(ct, rxn_id)
    ct.fluxes[rxn_id].bound = C.Between(-1000, 1000)
    ct.fluxes_forward[rxn_id].bound = C.Between(0, 1000)
    ct.fluxes_reverse[rxn_id].bound = C.Between(0, 1000)
    ct.interface.exchanges[rxn_id].bound = C.Between(-1000, 1000)
end

# fix some model quirks
function fix_bounds_ecfbc!(ct)
    ct.fluxes[:EX_5dglcn_e].bound = C.EqualTo(0.0)
    ct.fluxes[:EX_for_e].bound = C.EqualTo(0.0)
    ct.fluxes[:EX_pyr_e].bound = C.EqualTo(0.0)
end

# ##  Assemble community models

# load the basis model
model = load_model("iML1515.json", A.CanonicalModel.Model)

# identify membrane reactions (transporting stuff between compartments), used in the capacity bounds
membrane_rids = [
    rid for (rid, r) in model.reactions if
    length(unique(last.(split.(keys(r.stoichiometry), "_")))) != 1
]

# these co-cultures will be simulated
aa_pairs = [
    (:metA, :thrC),
    (:ilvA, :pheA),
    (:argA, :lysA),
    (:metA, :proA),
    (:argA, :metA),
    (:argA, :ilvA),
    (:metA, :pheA),
    (:ilvA, :lysA),
    (:ilvA, :metA),
]

# abundance range to screen through
abundances = [(a, 1 - a) for a in range(0.001, 0.999, length = 10)]

# parameters fed into screen function below
specs = [
    [(aa1, abun1), (aa2, abun2)] for (aa1, aa2) in aa_pairs for (abun1, abun2) in abundances
]

# ## Simulate cFBA

wt = flux_balance_constraints(model; interface = :identifier_prefixes)
open_bounds_fbc!(wt, :EX_glc__D_e)

cfba_res = screen(specs, workers = [1]) do spec
    auxotrophe_fba(wt, spec; fbc_only = true)
end

# ## Simulate ec-cFBA

wt = simplified_enzyme_constrained_flux_balance_constraints(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = [("membrane", membrane_rids, 166.0), ("total", A.reactions(model), 500.0)],
    interface = :identifier_prefixes,
)
open_bounds_ecfbc!(wt, :EX_glc__D_e)
fix_bounds_ecfbc!(wt)

eccfba_res = screen(specs. workers = [1]) do spec
    auxotrophe_fba(wt, spec)
end

#=

# ## Plot cFBA vs. ec-cFBA vs. experimental data
# Now that we have simulated the models, we need to compare them to data. The
# data from the paper is processed below, and then plotted.

using CairoMakie

# data from Mee et al., 2014
observed_abundances = [ # these pairings were simulated
    (:metA, :thrC, 17.0 / (17.0 + 80.8)),
    (:ilvA, :pheA, 15.2 / (15.2 + 51.0)),
    (:argA, :lysA, 11.7 / (11.7 + 37.1)),
    (:metA, :proA, 13.7 / (13.7 + 36.5)),
    (:argA, :metA, 22.8 / (22.8 + 19.8)),
    (:argA, :ilvA, 13.2 / (13.2 + 10.4)),
    (:metA, :pheA, 25.9 / (25.9 + 19.5)),
    (:ilvA, :lysA, 35.7 / (35.7 + 18.5)),
    (:ilvA, :metA, 71.3 / (71.3 + 17.3)),
]
observed_abundances = last.(observed_abundances) # in order of aa_pairs

# get maximum growth for each pairing in both simulations
cfba_abs = zeros(9) # optimal abundance using cFBA
eccfba_abs = zeros(9) # optimal abundance using ec-cFBA
for (i, (ko1, ko2)) in enumerate(aa_pairs)
    x = argmax(x -> x.mu, filter(x -> haskey(x, ko1) && haskey(x, ko2), cfba_res))
    cfba_abs[i] = x[ko1] / (x[ko1] + x[ko2])

    x = argmax(x -> x.mu, filter(x -> haskey(x, ko1) && haskey(x, ko2), eccfba_res))
    eccfba_abs[i] = x[ko1] / (x[ko1] + x[ko2])
end

# The figure below shows that using enzyme constraints dramatically improves the
# predictive validity of community simulations
fig = Figure();
ax = Axis(
    fig[1, 1],
    xlabel = "Observed composition",
    ylabel = "Predicted composition",
    title = "cFBA",
)
stem!(
    ax,
    observed_abundances,
    cfba_abs,
    offset = observed_abundances,
    marker = :rect,
    color = :blue,
    stemcolor = :blue,
    trunkcolor = :grey,
    label = "cFBA",
)
xlims!(ax, 0, 1)
ylims!(ax, 0, 1)

ax2 = Axis(fig[1, 2], xlabel = "Observed composition", title = "ec-cFBA")
stem!(
    ax2,
    observed_abundances,
    eccfba_abs,
    offset = observed_abundances,
    marker = :circle,
    color = :red,
    stemcolor = :red,
    trunkcolor = :grey,
    label = "ec-cFBA",
)
xlims!(ax2, 0, 1)
ylims!(ax2, 0, 1)
fig

=#
