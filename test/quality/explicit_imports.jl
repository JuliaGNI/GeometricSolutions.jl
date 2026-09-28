using ExplicitImports
using GeometricSolutions
using Test

test_explicit_imports(
    GeometricSolutions;
    # The package loads its dependencies with `using`, not with an explicit import list.
    no_implicit_imports = false,
    # Some imported names are not declared `public` by their owners.
    all_explicit_imports_are_public = false,
    # Some qualified accesses reach names that are not declared `public` by their owners.
    all_qualified_accesses_are_public = false
)
