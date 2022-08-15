push!(LOAD_PATH,"../src/")

using Documenter, MRIGradients

makedocs(sitename="MRIGradients Documentation",
    modules = [MRIGradients],
    pages = [
        "Home" => "index.md",
        "GIRFEssential" => "GIRFEssential.md",
        "GIRFApplier" => "GIRFApplier.md",
        "Utilities" => "Utilities.md"
    ]
)

deploydocs(;
    repo="github.com/BRAIN-TO/MRIGradients.jl",push_preview = true,
)
