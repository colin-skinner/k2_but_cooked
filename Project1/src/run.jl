using Project1: import_data, compute, write_gph, bayesian_score
using Graphs, GraphRecipes, Plots


#######################################################
#   Args
#######################################################

if length(ARGS) != 3
    error("usage: julia src/run.jl data/<infile>.csv [k2 retries] [local optimization runs]")
end

inputfilename = ARGS[1]
retries = parse(Int,ARGS[2])
local_iterations = parse(Int, ARGS[3])

cache_name = split(inputfilename, "/")[end]
cache_name = "cache/" * split(cache_name, ".")[1] * ".gph"

plot_filename = split(inputfilename, "/")[end]
plot_filename = "plots/" * split(plot_filename, ".")[1] # plot automatically adds png

outputfilename = split(inputfilename, "/")[end]
outputfilename = "graphs/" * split(outputfilename, ".")[1] * ".gph"

#######################################################
#   Importing 
#######################################################

nodes, samples = import_data(inputfilename)

old_graph = DiGraph() # default graph
order = Vector(range(1,length(nodes)))
old_score = -Inf

# Imports graph from cache if it exists
if isfile(cache_name)

    # Grab file and import header and data
    imported_graphs, order = open(cache_name, "r") do io
        order = eval(Meta.parse(readline(io)))
        graph = loadgraphs(io, LGFormat())
        return graph, order
    end

    if length(imported_graphs) == 0
        println("Imported file has no graphs. Defaulting to empty.")

    else
        old_graph = imported_graphs["graph"]
        old_score = bayesian_score(nodes, old_graph, samples)
        graph = old_graph
        
        println("Imported graph with $(nv(old_graph)) nodes")
        println("with score of $(bayesian_score(nodes, old_graph, samples))")
        println("and node order of $order")
    end
end

#######################################################
#   Importing 
#######################################################

# Computes new graph with score and order

graph, score, order = compute(nodes, samples, graph, retries, local_iterations)
println("Calculated score of $score")

# #######################################################
# #   Saving Cache
# #######################################################

# Save if new score is better

if score > old_score
    println("New graph score of ($score) is better than old score ($old_score). Saving...")
    open(cache_name, "w") do io
        println(io, order)
        savegraph(io, graph, "graph", LGFormat())
    end
else
    println("New graph score of ($score) was not better than old score ($old_score)")
end

#######################################################
#   Saving Graph 
#######################################################

node_names = [string(n.name) for n in nodes]

write_gph(graph, node_names, outputfilename)

graphplot(graph, names=node_names,size=(2000, 2000), dpi=300,
          nodesize=0.2,
          fontsize=12)
png(plot_filename)