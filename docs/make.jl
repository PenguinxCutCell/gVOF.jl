using Documenter
using gVOF

makedocs(
    modules = [gVOF],
    authors = "Fastaxx and contributors",
    repo = "https://github.com/PenguinxCutCell/gVOF.jl/blob/{commit}{path}#{line}",
    sitename = "gVOF.jl",
    format = Documenter.HTML(
        canonical = "https://PenguinxCutCell.github.io/gVOF.jl",
        repolink = "https://github.com/PenguinxCutCell/gVOF.jl",
        collapselevel = 2,
    ),
    pages = [
        "Home" => "index.md",
        "Types" => "types.md",
        "Grid" => "grid.md",
        "Tagging and Initialization" => "tag-init.md",
        "Reconstruction" => "reconstruction.md",
        "Advection" => "advection.md",
        "Output and Test Cases" => "output-testcases.md",
        "Reference" => "95-reference.md",
    ],
    pagesonly = true,
    warnonly = true,
)

deploydocs(
    repo = "github.com/PenguinxCutCell/gVOF.jl",
    push_preview = true,
)
