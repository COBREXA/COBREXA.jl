FROM julia

RUN julia -e 'import Pkg; Pkg.add(["COBREXA", "HiGHS", "Clarabel", "Tulip", "GLPK", "ConstraintTrees", "AbstractFBCModels", "JSONFBCModels", "MATFBCModels", "SBMLFBCModels"]); Pkg.resolve(); Pkg.status(); Pkg.instantiate(); Pkg.precompile()'

CMD ["julia"]
