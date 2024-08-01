
<div align="center">
    <img src="docs/src/assets/header.svg?maxAge=0" width="80%">
</div>

# COnstraint-Based Reconstruction and EXascale Analysis

[docs-img-stable]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-url-stable]: https://cobrexa.github.io/COBREXA.jl

[docs-img-dev]: https://img.shields.io/badge/docs-latest-0af.svg
[docs-url-dev]: https://cobrexa.github.io/COBREXA.jl/dev/
[docs-url-examples]: https://cobrexa.github.io/COBREXA.jl/dev/examples/
[docs-url-fba]: https://cobrexa.github.io/COBREXA.jl/dev/examples/02a-flux-balance-analysis/

[docker-url]: https://github.com/COBREXA/COBREXA.jl/pkgs/container/cobrexa.jl
[docker-img]: https://ghcr-badge.egpl.dev/cobrexa/cobrexa.jl/size?color=%2362a0ea&tag=latest&label=docker&trim=

[ci-img]: https://github.com/COBREXA/COBREXA.jl/actions/workflows/ci.yml/badge.svg?branch=master
[ci-url]: https://github.com/COBREXA/COBREXA.jl/actions/workflows/ci.yml

[cov-img]: https://codecov.io/gh/COBREXA/COBREXA.jl/branch/master/graph/badge.svg?token=H3WSWOBD7L
[cov-url]: https://codecov.io/gh/COBREXA/COBREXA.jl

[legacy-url]: https://github.com/LCSB-BioCore/COBREXA.jl

[repostatus-url]: https://www.repostatus.org/#active
[repostatus-img]: https://www.repostatus.org/badges/latest/active.svg

| **Documentation** | **Tests** | **Coverage** | **Project status** |
|:---:|:---:|:---:|:---:|
| [![docs-img-stable]][docs-url-stable] [![docs-img-dev]][docs-url-dev] | [![CI][ci-img]][ci-url] | [![codecov][cov-img]][cov-url] | [![repostatus-img]][repostatus-url] |

COBREXA.jl provides constraint-based reconstruction and analysis tools for
exa-scale metabolic modeling in [Julia](https://julialang.org/).

For COBREXA.jl releases before version 2.0 (released in 2023 and earlier),
please see [the legacy COBREXA repository with version 0.x and 1.x][legacy-url].

## Why use COBREXA.jl?

Compared to the other modeling toolkits, COBREXA.jl offers:

- Almost effortless parallelization to scale-up to HPC environments. This
  allows the users to quickly run huge analyses.
- A straightforward way to construct new kinds of constraint systems by
  re-using and re-combining the existing ones with newly implemented
  functionality, provided by
  [ConstraintTrees.jl](https://github.com/COBREXA/ConstraintTrees.jl)
  middleware. This vastly simplifies reconstruction of complex systems, such as
  multi-compartment ("community") growth-balanced models with both per-member
  and global resource-allocation and metabolite availability constraints.
- Performance at the level of "raw" user workflow operations, because all
  analysis code are compiled and optimized by Julia. This reduces the
  computation time and memory required to analyze very large models.
- A wide range of pre-implemented HPC-compatible analysis functions, with
  [tested examples available in the documentation][docs-url-examples].

## Getting started

[COBREXA.jl documentation][docs-url-stable]
is available online (you may also be interested in the
[development version][docs-url-dev]
of the package documentation).

### Installation

COBREXA.jl is best installed using Julia's packaging system.

1. If you do not have Julia installed, [download it from the official
   site](https://julialang.org/downloads/).
2. Start Julia by typing `julia` into a terminal
3. Type `] add COBREXA`. This will install COBREXA automatically, with all
   dependencies.
4. You may want to install [additional constraint-solver software compatible
   with JuMP](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers).
   Type e.g. `] add GLPK` to install the popular GLPK solver.

### Guides

To start quickly, [follow the examples in the documentation][docs-url-examples].

We recommend to start with the usual [Flux Balance Analysis example][docs-url-fba].

### Testing the installation

If you run a non-standard platform (e.g. a customized operating system), or if
you added any modifications to the `COBREXA` source code, you may want to run
the test suite to ensure that everything works as expected:

```julia
] test COBREXA
```

### Prebuilt Docker images [![docker][docker-img]][docker-url]

Docker image is available from the docker hub as
[cobrexa/cobrexa.jl][docker-url], and from GitHub container repository.
Download and use them as usual with docker:

```sh
docker run -ti --rm cobrexa/cobrexa.jl:latest

# or alternatively from ghcr.io
docker run -ti --rm ghcr.io/cobrexa/docker/cobrexa.jl:latest
```

In the container, you should get a `julia` shell with the important packages
already installed, and you may immediately import the package with `using
COBREXA`, and start running COBREXA analyses (such as the ones in the
[documentation examples][docs-url-examples]).

If you require precise reproducibility, use a tag like `v1.2.2` instead of
`latest` (all releases since 1.2.2 are tagged this way).

**NOTE:** Versions of COBREXA before 2.0 are available from the original
repository at [github.com/LCSB-BioCore/COBREXA.jl][legacy-url].

## Acknowledgements

`COBREXA.jl` is developed at the Luxembourg Centre for Systems Biomedicine of
the University of Luxembourg ([uni.lu/lcsb](https://wwwen.uni.lu/lcsb)),
cooperating with the Institute for Quantitative and Theoretical Biology at the Heinrich
Heine University in Düsseldorf ([qtb.hhu.de](https://www.qtb.hhu.de/)).

The development was supported by European Union's Horizon 2020 Programme under
PerMedCoE project ([permedcoe.eu](https://permedcoe.eu/)) agreement no. 951773.

If you use COBREXA.jl and want to refer to it in your work, use the following
citation format (also available as BibTeX in [cobrexa.bib](cobrexa.bib)):

> Miroslav Kratochvíl, Laurent Heirendt, St Elmo Wilken, Taneli Pusa, Sylvain Arreckx, Alberto Noronha, Marvin van Aalst, Venkata P Satagopam, Oliver Ebenhöh, Reinhard Schneider, Christophe Trefois, Wei Gu, *COBREXA.jl: constraint-based reconstruction and exascale analysis*, Bioinformatics, Volume 38, Issue 4, 15 February 2022, Pages 1171–1172, https://doi.org/10.1093/bioinformatics/btab782

<img src="docs/src/assets/cobrexa.svg" alt="COBREXA logo" height="64px" style="height:64px; width:auto">   <img src="docs/src/assets/unilu.svg" alt="Uni.lu logo" height="64px" style="height:64px; width:auto">   <img src="docs/src/assets/lcsb.svg" alt="LCSB logo" height="64px" style="height:64px; width:auto">   <img src="docs/src/assets/hhu.svg" alt="HHU logo" height="64px" style="height:64px; width:auto">   <img src="docs/src/assets/qtb.svg" alt="QTB logo" height="64px" style="height:64px; width:auto">   <img src="docs/src/assets/permedcoe.svg" alt="PerMedCoE logo" height="64px" style="height:64px; width:auto">
