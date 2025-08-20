# ! https://github.com/MakieOrg/Makie.jl/blob/master/docs/attrdocs_block.jl

abstract type AttrdocsBlocks <: Documenter.Expanders.NestedExpanderPipeline end

Documenter.Selectors.order(::Type{AttrdocsBlocks}) = 8.0 # like @example
function Documenter.Selectors.matcher(::Type{AttrdocsBlocks}, node, page, doc)
    Documenter.iscode(node, r"^@attrdocs")
end

# this type is just for a helper node which we can push child elements to that are
# at the end flattened as if they were all on the level of the container
struct Container <: Documenter.AbstractDocumenterBlock
    codeblock::MarkdownAST.CodeBlock
end

MarkdownAST.can_contain(::Container, ::MarkdownAST.AbstractElement) = true

function DocumenterVitepress.render(io::IO, mime::MIME"text/plain", node::MarkdownAST.Node,
                                    c::Container, page, doc; kwargs...)
    return DocumenterVitepress.render(io, mime, node, node.children, page, doc; kwargs...)
end

function attrs_examples_docs_defaults(type::Type{<:Makie.Block})
    attrkeys = sort(collect(keys(Makie.default_attribute_values(type, nothing))))

    all_examples = Makie.attribute_examples(type)
    all_docs = Makie._attribute_docs(type)
    all_defaults = Makie.attribute_default_expressions(type)

    return (; attrkeys, all_examples, all_docs, all_defaults)
end

function attrs_examples_docs_defaults(type::Type{<:Makie.Plot})
    docatt = Makie.documented_attributes(type)
    metadata = docatt.d

    attrkeys = sort(collect(keys(metadata)))
    all_examples = Makie.attribute_examples(type)
    all_docs = Dict([(attr => something(meta.docstring, "No docs available."))
                     for (attr, meta) in metadata])
    all_defaults = Dict([(attr => meta.default_expr) for (attr, meta) in metadata])

    return (; attrkeys, all_examples, all_docs, all_defaults)
end

function Documenter.Selectors.runner(::Type{AttrdocsBlocks}, node, page, doc)
    codeblock = node.element
    node.element = Container(codeblock)

    type = getproperty(TimeseriesTools, Symbol(strip(codeblock.code)))

    (; attrkeys, all_examples, all_docs, all_defaults) = attrs_examples_docs_defaults(type)

    for attrkey in attrkeys
        heading = @ast MarkdownAST.Heading(3) do
            "$attrkey"
        end
        push!(node.children, heading)

        default_str = all_defaults[attrkey]
        default_para = @ast MarkdownAST.Paragraph() do
            "Defaults to "
            MarkdownAST.Code(default_str)
        end
        push!(node.children, default_para)

        docs = get(all_docs, attrkey, nothing)
        docs_md = Markdown.parse(docs)
        docs_node = convert(MarkdownAST.Node, docs_md)
        append!(node.children, docs_node.children)

        examples = get(all_examples, attrkey, Makie.Example[])

        for example in examples
            figurenode = @ast MarkdownAST.CodeBlock("@figure", example.code)
            push!(node.children, figurenode)
        end
    end

    # expand the children we added with the normal expand pipeline, for example to apply the standard
    # treatment to headings that documenter would normally do (which does not happen in the nested expand pipeline)
    # and for running the @figure nodes
    for childnode in collect(node.children)
        Documenter.Selectors.dispatch(Documenter.Expanders.ExpanderPipeline, childnode,
                                      page, doc)
        Documenter.expand_recursively(childnode, page, doc)
    end

    return
end

# ! https://github.com/MakieOrg/Makie.jl/blob/master/docs/shortdocs_block.jl

abstract type ShortDocsBlocks <: Documenter.Expanders.NestedExpanderPipeline end

Documenter.Selectors.order(::Type{ShortDocsBlocks}) = 3.0 # like @docs
function Documenter.Selectors.matcher(::Type{ShortDocsBlocks}, node, page, doc)
    Documenter.iscode(node, r"^@shortdocs")
end

function unlink_with_all_following_siblings!(node)
    next = node.next
    MarkdownAST.unlink!(node)
    if next !== nothing
        unlink_with_all_following_siblings!(next)
    end
    return
end

# ```@shortdocs is like ```@docs but it cuts off the docstring for the plot functions at the attributes
# section because all the attributes are already added with examples anyway
function Documenter.Selectors.runner(::Type{ShortDocsBlocks}, node, page, doc)
    el = node.element
    el.info = replace(el.info, "@shortdocs" => "@docs")
    Documenter.Selectors.runner(Documenter.Expanders.DocsBlocks, node, page, doc)

    docsnode = first(node.children).element
    if !(docsnode isa Documenter.DocsNode)
        error("docs node conversion failed for $el")
    end

    mdasts = docsnode.mdasts

    ast_to_look_for = MarkdownAST.@ast MarkdownAST.Paragraph() do
        MarkdownAST.Strong() do
            "Attributes"
        end
    end

    for mdast in mdasts
        found = false
        for child in mdast.children
            if child == ast_to_look_for
                unlink_with_all_following_siblings!(child)
                found = true
                break
            end
        end
        if !found
            display(mdast)
            error("Found no Attributes section in above markdown ast")
        end
    end
    return
end
