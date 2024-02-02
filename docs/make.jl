using TimeseriesTools
using Unitful
using Documenter

DocMeta.setdocmeta!(TimeseriesTools, :DocTestSetup, :(using Unitful, TimeseriesTools);
                    recursive = true)

makedocs(;
         modules = [TimeseriesTools],
         authors = "brendanjohnharris <brendanjohnharris@gmail.com> and contributors",
         repo = "https://github.com/brendanjohnharris/TimeseriesTools.jl/blob/{commit}{path}#{line}",
         sitename = "TimeseriesTools.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://brendanjohnharris.github.io/TimeseriesTools.jl",
                                  edit_link = "main",
                                  assets = String[],),
         pages = [
             "Home" => "index.md",
             "Types" => "types.md",
             "Utils" => "utils.md",
             "Others" => "others.md",
         ],)

deploydocs(;
           repo = "github.com/brendanjohnharris/TimeseriesTools.jl",
           devbranch = "main",)
