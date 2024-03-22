
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

Helper for getting stuff from dictionaries where keys may be easily missing.
"""
maybeget(_::Nothing, _...) = nothing
maybeget(x, k, ks...) = haskey(x, k) ? maybeget(x[k], ks...) : nothing
maybeget(x) = x

"""
$(TYPEDSIGNATURES)

Helper for applying functions to stuff that might be `nothing`.
"""
maybemap(f, ::Nothing, def = nothing) = def
maybemap(f, x, ::Any = nothing) = f(x)
