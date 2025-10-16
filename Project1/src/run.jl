using Project1: import_data, compute, write_gph, bayesian_score
using Graphs


#######################################################
#   Importing 
#######################################################

if length(ARGS) != 1
    error("usage: julia src/run.jl data/<infile>.csv")
end

inputfilename = ARGS[1]
# outputfilename = ARGS[2]
cache_name = split(inputfilename, "/")[end]
cache_name = "cache/" * split(cache_name, ".")[1] * ".gph"

outputfilename = split(inputfilename, "/")[end]
outputfilename = "graphs/" * split(outputfilename, ".")[1] * ".gph"

nodes, samples = import_data(inputfilename, false)

old_graph = DiGraph() # default graph

# Imports graph from cache if it exists
if isfile(cache_name)
    imported_graphs = open(cache_name, "r") do io
        loadgraphs(io, LGFormat())
    end

    if length(imported_graphs) == 0
        println("Imported file has no graphs. Defaulting to empty.")
    else
        old_graph = imported_graphs["graph"]
        println("Imported graph with $(nv(old_graph)) nodes")
        println("with score of $(bayesian_score(nodes, old_graph, samples))")
    end
end

println("Number of old_graph nodes: $(nv(old_graph))")

#######################################################
#   Importing 
#######################################################

# Computes new graph with score and order
new_graph, score, sorted_nodes, new_nodes, new_samples = compute(nodes, samples, old_graph)

println("Calculated score of $score")

#######################################################
#   Saving Cache
#######################################################

# Saves graph
if nv(old_graph) == 0
    println("Saving graph with score of $score")
    open(cache_name, "w") do io
        savegraph(io, new_graph, "graph", LGFormat())
    end
elseif score > bayesian_score(new_nodes, old_graph, new_samples)
    println("New graph score of ($score) is better than old score ($(bayesian_score(new_nodes, old_graph, new_samples))). Saving...")
    open(cache_name, "w") do io
        savegraph(io, new_graph, "graph", LGFormat())
    end
else
    println("New graph score of ($score) was not better than old score ($(bayesian_score(new_nodes, old_graph, new_samples)))")

end

#######################################################
#   Saving Graph 
#######################################################

# New order of names
idx2names = Matrix{Symbol}(undef, 1, length(nodes))
for (node, new_idx) in zip(nodes, sorted_nodes.ordering)
    idx2names[new_idx] = node.name
end

write_gph(new_graph, idx2names, outputfilename)


