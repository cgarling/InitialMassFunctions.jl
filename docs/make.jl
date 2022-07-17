using Documenter, Distributions, InitialMassFunctions

# The `format` below makes it so that urls are set to "pretty" if you are pushing them to a hosting service, and basic if you are just using them locally to make browsing easier.

makedocs(
    sitename="InitialMassFunctions.jl",
    modules = [InitialMassFunctions],
    format = Documenter.HTML(;prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Chris Garling",
    pages = ["index.md","api.md"],
    doctest=true
)
