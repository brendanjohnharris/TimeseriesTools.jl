import FileIO: save, load, @format_str, File, query
using JSON
using DelimitedFiles

export savetimeseries, savets, loadtimeseries, loadts

savetimeseries(f::String, x) = savetimeseries(f|>query, x)
loadtimeseries(f::String) = loadtimeseries(f|>query)

## JLD2 files are easiest
savetimeseries(f::File{format"JLD2"}, x::AbstractTimeSeries) = save(f, Dict("timeseries" => x ))
loadtimeseries(f::File{format"JLD2"}) = load(f, "timeseries")

## CSV files are harder. We can't fully reconstruct a generic timeseries, so need to make some concessions.
# We'll assume that the first column is the time index, and the remaining columns are the data.
function savetimeseries(f::File{format"CSV"}, x::AbstractTimeSeries, var)
    isnothing(var) && (var = "")
    open(f.filename, "w") do f
        print(f, "# ")

        if name(x) isa DimensionalData.NoName
            print(f, "")
        else
            try
                print(f, json(name(x)))
            catch e
                @warn e
            end
        end

        print(f, "\n# ")

        if metadata(x) isa DimensionalData.Dimensions.LookupArrays.NoMetadata
            print(f, "")
        else
            try
                print(f, json(metadata(x)))
            catch e
                @warn e
            end
        end

        print(f, "\n# ")
        try
            refs = isempty(refdims(x)) ? "" : json(refdims(x))
            print(f, refs)
        catch e
            print(f, "")
            @warn e
        end

        print(f, "\n# ")

        try
            vars = isnothing(dims(x, Var)) ? ["time"] : json(["time", var.val.data])
        catch e
            vars = json(["time", 1:length(var)]) # egh
            @warn e
        end


        print(f, vars)
        vars = join(var, ',')
        print(f, "\ntime,$vars\n")
        writedlm(f, [times(x) x.data], ',')
    end
end
function savetimeseries(f::File{format"CSV"}, x::UnivariateTimeSeries)
    if length(refdims(x)) == 1
        var = refdims(x)
    else
        var = refdims(x, Var)
    end
    savetimeseries(f, x, var)
end
savetimeseries(f::File{format"CSV"}, x::MultivariateTimeSeries) = savetimeseries(f, x, dims(x, 2))

function loadtimeseries(f::File{format"CSV"})
    open(f.filename, "r") do f
        # Read the name
        line = readline(f)
        name = isempty(line[3:end]) ? DimensionalData.NoName() : JSON.parse(line[3:end])

        # Read the metadata
        line = readline(f)
        metadata = isempty(line[3:end]) ? DimensionalData.Dimensions.LookupArrays.NoMetadata() : JSON.parse(line[3:end])

        # Read the reference dimensions
        line = readline(f)
        if !isempty(line[3:end])
            @warn "Cannot load refdims yet"
        end
        refdims = () # Not implemented, maybe isempty(line[3:end]) ? () : JSON.parse(line[3:end])

        # Read the variables
        line = readline(f)
        j = JSON.parse(line[3:end])
        vars = length(j) > 1 ? j[2] : ()

        data, _ = readdlm(f, ',', header=true)

        if isempty(vars)
            x = TimeSeries(data[:, 1], data[:, 2]; name, metadata, refdims)
        else
            x = TimeSeries(data[:, 1], vars, data[:, 2:end]; name, metadata, refdims)
        end
    end
end

savets = savetimeseries
loadts = loadtimeseries
