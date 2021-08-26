using IndirectCoStateInit
using Documenter

DocMeta.setdocmeta!(IndirectCoStateInit, :DocTestSetup, :(using IndirectCoStateInit); recursive=true)

makedocs(;
    modules=[IndirectCoStateInit],
    authors="Grant Hecht",
    repo="https://github.com/GrantHecht/IndirectCoStateInit.jl/blob/{commit}{path}#{line}",
    sitename="IndirectCoStateInit.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://GrantHecht.github.io/IndirectCoStateInit.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/GrantHecht/IndirectCoStateInit.jl",
    devbranch="main",
)
