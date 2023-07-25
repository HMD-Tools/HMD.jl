using HMD
using Documenter

DocMeta.setdocmeta!(HMD, :DocTestSetup, :(using HMD); recursive=true)

makedocs(;
    modules=[HMD],
    authors="yoshiatsu163 <atsushi.yoshida.github@gmail.com> and contributors",
    repo="https://github.com/yoshiatsu163/HMD.jl/blob/{commit}{path}#{line}",
    sitename="HMD.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yoshiatsu163.github.io/HMD.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/yoshiatsu163/HMD.jl",
    devbranch="master",
)
