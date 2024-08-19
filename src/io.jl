
# Copyright (c) 2021-2024, University of Luxembourg
# Copyright (c) 2021-2024, Heinrich-Heine University Duesseldorf
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
$(TYPEDSIGNATURES)

Load a FBC model representation while guessing the correct model type to load.
Uses `AbstractFBCModels.load`.

This overload almost always involves a search over types; do not use it in
environments where performance is critical.
"""
function load_model(path::String)
    A.load(path)
end

"""
$(TYPEDSIGNATURES)

Overload of [`load_model`](@ref) that guesses the input type, but immediately
converts to the model type given by argument `convert_to`.
"""
function load_model(path::String, convert_to::Type{O}) where {O<:A.AbstractFBCModel}
    convert(O, A.load(path))
end

"""
$(TYPEDSIGNATURES)

Load a FBC model representation from a known `model_type`. Uses
`AbstractFBCModels.load`.
"""
function load_model(model_type::Type{I}, path::String) where {I<:A.AbstractFBCModel}
    A.load(I, path)
end

"""
$(TYPEDSIGNATURES)

Overload of [`load_model`](@ref) that explicitly specifies the known input
type, and immediately converts to another model type given by argument
`convert_to`.
"""
function load_model(
    model_type::Type{I},
    path::String,
    convert_to::Type{O},
) where {I<:A.AbstractFBCModel,O<:A.AbstractFBCModel}
    convert(O, A.load(I, path))
end

export load_model

"""
$(TYPEDSIGNATURES)

Save a FBC model representation. Uses `AbstractFBCModels.save`.

Use the 3-parameter overload if you need to convert the model to another
representation (e.g., if you want to save a canonical model type as JSON or
SBML).
"""
function save_model(model::T, path::String) where {T<:A.AbstractFBCModel}
    A.save(model, path)
end

"""
$(TYPEDSIGNATURES)

Overload of [`save_model`](@ref) that converts the model type to `convert_to`
before saving.
"""
function save_model(
    model::T,
    path::String,
    convert_to::Type{O},
) where {T<:A.AbstractFBCModel,O<:A.AbstractFBCModel}
    A.save(convert(O, model), path)
end

export save_model

"""
$(TYPEDSIGNATURES)

Like [`save_model`](@ref) but tries to convert the `model` to a type that
matches the extension of the `path`. For example, this will convert the
`model` to a JSON model type in case the `path` ends with `.json`.

This is an utility shortcut -- if possible, it is always better to specify the
output model type explicitly.
"""
function save_converted_model(model::T, path::String) where {T<:A.AbstractFBCModel}
    save_model(model, path, A.guess_model_type_from_filename(path))
end

export save_converted_model

"""
$(TYPEDSIGNATURES)

Safely download a model with a known hash. All arguments are forwarded to
`AbstractFBCModels.download_data_file` -- see the documentation in the
AbstractFBCModels package for details.
"""
download_model(args...; kwargs...) = A.download_data_file(args...; kwargs...)

export download_model
