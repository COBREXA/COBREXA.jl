
# Copyright (c) 2021-2025, University of Luxembourg
# Copyright (c) 2021-2025, Heinrich-Heine University Duesseldorf
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
    module COBREXA

COnstraint Based Reconstruction and EXascale Analysis. COBREXA provides
functions for construction, modification, simulation and analysis of
constraint-based metabolic models that follows the COBRA methodology.

COBREXA is built as a front-end for the combination of `AbstractFBCModels.jl`
(provides the model I/O), `ConstraintTrees.jl` (provides the constraint system
organization), `Distributed.jl` (provides HPC execution capability), and
`JuMP.jl` (provides the solvers).

See the online documentation for a complete description of functionality aided
by copy-pastable examples.

To start quickly, load your favorite JuMP-compatible solver, use
[`load_model`](@ref) to read a metabolic model from the disk, and solve it with
[`flux_balance_analysis`](@ref).
"""
module COBREXA

using DocStringExtensions

import AbstractFBCModels as A
import ConstraintTrees as C
import Distributed as D
import JuMP as J
import LinearAlgebra
import SparseArrays
import StableRNGs

include("types.jl")

# core functionality
include("io.jl")
include("solver.jl")
include("worker_data.jl")

# generic analysis functions
include("analysis/envelope.jl")
include("analysis/frontend.jl")
include("analysis/parsimonious.jl")
include("analysis/sample.jl")
include("analysis/screen.jl")
include("analysis/solver.jl")
include("analysis/threshold.jl")
include("analysis/variability.jl")

# conversion of various stuff to constraint trees
include("builders/compare.jl")
include("builders/enzymes.jl")
include("builders/interface.jl")
include("builders/knockout.jl")
include("builders/loopless.jl")
include("builders/objectives.jl")
include("builders/scale.jl")
include("builders/unsigned.jl")

# simplified front-ends for the above
include("frontend/balance.jl")
include("frontend/community.jl")
include("frontend/concentrations.jl")
include("frontend/envelope.jl")
include("frontend/enzymes.jl")
include("frontend/gapfill.jl")
include("frontend/knockout.jl")
include("frontend/loopless.jl")
include("frontend/medium.jl")
include("frontend/mmdf.jl")
include("frontend/moma.jl")
include("frontend/parsimonious.jl")
include("frontend/sample.jl")
include("frontend/variability.jl")

# utilities
include("misc/bounds.jl")
include("misc/breaks.jl")
include("misc/maybe.jl")
include("misc/settings.jl")
include("misc/trees.jl")

# configuration (initializes with values from other source files, thus last)
include("config.jl")

end # module COBREXA
