
"""
    abstract type AbstractMetabolicModel end

A helper supertype of everything usable as a linear-like model for COBREXA
functions.

If you want your model type to work with COBREXA, add the `AbstractMetabolicModel` as
its supertype, and implement the accessor functions. Accessors
[`reactions`](@ref), [`metabolites`](@ref), [`stoichiometry`](@ref),
[`bounds`](@ref) and [`objective`](@ref) must be implemented; others are not
mandatory and default to safe "empty" values.
"""
abstract type AbstractMetabolicModel end

"""
    abstract type AbstractModelWrapper <: AbstractMetabolicModel end

A helper supertype of all "wrapper" types that contain precisely one other
[`AbstractMetabolicModel`](@ref).
"""
abstract type AbstractModelWrapper <: AbstractMetabolicModel end

const SparseMat = SparseMatrixCSC{Float64,Int}
const SparseVec = SparseVector{Float64,Int}
const MatType = AbstractMatrix{Float64}
const VecType = AbstractVector{Float64}
const StringVecType = AbstractVector{String}

"""
    MetaboliteFormula = Dict{String,Int}

Dictionary of atoms and their abundances in a molecule.
"""
const MetaboliteFormula = Dict{String,Int}

"""
    Annotations = Dict{String,Vector{String}}

Dictionary used to store (possible multiple) standardized annotations of
something, such as a [`Metabolite`](@ref) and a [`Reaction`](@ref).

# Example
```
Annotations("PubChem" => ["CID12345", "CID54321"])
```
"""
const Annotations = Dict{String,Vector{String}}

"""
    Notes = Dict{String,Vector{String}}

Free-form notes about something (e.g. a [`Gene`](@ref)), categorized by
"topic".
"""
const Notes = Dict{String,Vector{String}}

"""
    GeneAssociationsDNF = Vector{Vector{String}}

Disjunctive normal form of simple gene associations. For example, `[[A, B],
[B]]` represents two isozymes where the first requires both genes `A` and `B`,
while the second isozyme only requires gene `C`.

This string representation is typically used to represent gene reaction rules,
but does not contain any subunit stoichiometry of kinetic information of the
isozymes. See [`Isozyme`}(@ref) for a more complete structure.
"""
const GeneAssociationsDNF = Vector{Vector{String}}