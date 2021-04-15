using Documenter, GuaranteedEllipseFitting

format = Documenter.HTML(edit_link = "master",
                         prettyurls = get(ENV, "CI", nothing) == "true",
                         assets = String[])

makedocs(;
    modules=[GuaranteedEllipseFitting],
    format = format,
    pages=[
        "Home" => "index.md",
        "Function Reference" => "reference.md"
    ],
    sitename="GuaranteedEllipseFitting.jl",
    authors="Dr. Zygmunt L. Szpak"
)

deploydocs(;
    repo="github.com/JuliaImages/GuaranteedEllipseFitting.jl.git",
)
