using CairoMakie
using Makie
import Makie.Linestyle
using TimeseriesTools
using TimeseriesMakie
using Documenter
using Documenter: Documenter
using Documenter.MarkdownAST
using Documenter.MarkdownAST: @ast
using DocumenterVitepress
using Markdown

ENV["UNITFUL_FANCY_EXPONENTS"] = true

include("docs_blocks.jl")

format = DocumenterVitepress.MarkdownVitepress(;
                                               repo = "github.com/brendanjohnharris/TimeseriesTools.jl",
                                               devbranch = "main",
                                               devurl = "dev")

timeseriestools = [
    "Introduction" => "TimeseriesTools/index.md",
    "Types" => "TimeseriesTools/types.md",
    "Utils" => "TimeseriesTools/utils.md",
    "Others" => "TimeseriesTools/others.md"
]

timeseriesmakie = ["Introduction" => "TimeseriesMakie/index.md",
    "Recipes" => "TimeseriesMakie/recipes.md",
    "Reference" => "TimeseriesMakie/reference.md"]

pages = ["Home" => "index.md",
    "Quick start" => "quickstart.md",
    "TimeseriesTools" => timeseriestools,
    "TimeseriesMakie" => timeseriesmakie
]

makedocs(;
         authors = "brendanjohnharris <brendanjohnharris@gmail.com> and contributors",
         sitename = "TimeseriesTools",
         format,
         doctest = false,
         warnonly = [:cross_references],
         modules = [TimeseriesTools, TimeseriesMakie],
         pages)

DocumenterVitepress.deploydocs(;
                               repo = "github.com/brendanjohnharris/TimeseriesTools.jl",
                               target = "build", # this is where Vitepress stores its output
                               branch = "gh-pages",
                               devbranch = "main",
                               push_preview = true)
