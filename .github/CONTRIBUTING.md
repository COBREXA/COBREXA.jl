# Contributing to COBREXA.jl

The package follows the usual GitHub workflow:
- Issues and discussions are welcome. Please try to follow the issue templates.
- Pull requests are welcome.

## Development hints

- use Julia's package-development mode (`] dev COBREXA` will give you a source
  tree to work with, or equivalently `] dev --local some/path` will load the package
  from some given path where you git-cloned the repository before)
- try to follow the package coding standards -- these include mainly
  - code functionally, avoid functions with side effects
  - code less (shorter code is always better)
  - make sure your functions compose well with the other ones in the package
  - keep full, unicode-free names of all functions and variables
  - all use-case-specific code should have a generalized version available to users
  - use JuliaFormatter
  - write docstrings, use DocStringExtensions
  - do not break existing tests (in Julia: `] test COBREXA`), unless you have a
    proof that breaking them is the correct thing to do
  - if possible, prefix branch names with your initials (`jd-something` is a
    good branch from some random John Doe)

## Mergeability conditions

We will welcome any contributions, but we strive to keep COBREXA highly
reliable and efficient. This limits the scope of the pull requets that we can
merge.

For the above reason, we might have to reject PRs that contain some of the
following:

- large new dependencies (especially ones that increase the package load time)
- bundled code or redundant functionality from other packages
- new features that are not covered with tested documentation (see
  `docs/src/examples`), and generally code that is not covered by tests
- very use-case-specific features that do not generalize very well (e.g.,
  model-specific analysis functions)
- code with unclear licensing terms or incompatible with the repository license
- removal of documented features

If you are not sure whether your code can be merged, or if there are problems
with any of the above points, please still feel free to open a PR, and leave a
notice -- we will see what we can do. Typically, all contributions will
eventually find a spot where they will be useful.

As a work-around for the above issues, it is often useful to consider creating
a new Julia package that extends COBREXA.
