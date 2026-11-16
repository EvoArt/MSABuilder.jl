using MSABuilder
using Documenter

makedocs(;
    modules=[MSABuilder],
    authors="Arthur Newbury",
    sitename="MSABuilder.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EvoArt.github.io/MSABuilder.jl",
        edit_link="master",
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
    ],
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/EvoArt/MSABuilder.jl",
    devbranch="master",
)
