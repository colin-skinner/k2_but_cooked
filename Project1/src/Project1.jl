module Project1

using Graphs
using Printf
using LinearAlgebra
using SpecialFunctions
using Random

#######################################################
#   Datatypes
#######################################################

mutable struct Variable
    name::Symbol
    r::Int # number of possible values
end

# TODO: decide if using
struct K2_Ordering
    ordering::Vector{Int} # variable ordering
end

struct LocalDirectedGraphSearch
    G::SimpleDiGraph # initial graph
    k_max::Int # number of iterations
end

#######################################################
#   Importing and Exporting
#######################################################

function import_data(filename::String, debug::Bool = false)::Tuple{Vector{Variable}, Matrix{Int}}

    nodes = Vector{Variable}()
    samples = Vector{Vector{Int}}()

    open(filename, "r") do io

        # Break if the file is empty
        header = readline(io)
        if eof(io)
            return nodes, Matrix{Int}(undef, 0, 0)
        end


        # Parse node names, setting default to 1 possible value
        node_names = split(header, ',')
        nodes = [Variable(Symbol(name), 1) for name in node_names]

        # Parse samples
        for line in eachline(io)
            single_sample = parse.(Int, split(line, ','))
            push!(samples, single_sample)
            
            # if debug
            #     @printf("%s\n", single_sample)
            # end
        end
    end

    # Convert samples to a Matrix
    if !isempty(samples)
        num_samples = length(samples)
        sample_matrix = Matrix{Int}(undef, length(nodes), num_samples)
        
        for i in 1:num_samples
            sample_matrix[:, i] = samples[i]
        end
    else
        sample_matrix = Matrix{Int}(undef, 0, 0)
    end

    # Max values
    for (i, var_values) in enumerate(eachrow(sample_matrix))
        nodes[i].r = maximum(var_values)
    end

    return nodes, sample_matrix
end

"""
    write_gph(dag::DiGraph, idx2names, filename)

Takes a DiGraph, a Dict of index to names and a output filename to write the graph in `gph` format.
"""
function write_gph(dag::DiGraph, idx2names, filename)
    open(filename, "w") do io
        for edge in edges(dag)
            @printf(io, "%s,%s\n", idx2names[src(edge)], idx2names[dst(edge)])
        end
    end
end

#######################################################
#   Utility Functions
#######################################################

# From book (with TA recommendation)
function sub2ind(siz, x)
    k = vcat(1, cumprod(siz[1:end-1]))
    return dot(k, x .- 1) + 1
end

function statistics(nodes::Vector{Variable}, graph::SimpleDiGraph, data::Matrix{Int})
    n = size(data, 1)
    r_i = [nodes[i].r for i in 1:n]

    q_i = [ # Parental instantiations for each node
        prod([r_i[i] for i in inneighbors(graph,node)]) 
        for node in 1:n]

    # List of counts matrices
    M = [zeros(q_i[i], r_i[i]) for i in 1:n]

    # Iterate through all each variable
    for node_data in eachcol(data)
        for node_num_i in 1:n
            k = node_data[node_num_i]
            parents = inneighbors(graph, node_num_i)

            parental_instantiations_j = 1       # Default case of no parents
            if length(parents) > 0
                parental_instantiations_j = sub2ind(r_i[parents], node_data[parents])
            end
            
            M[node_num_i][parental_instantiations_j, k] += 1
        end
    end

    return M
end

function prior(nodes::Vector{Variable}, graph::SimpleDiGraph)
    n = length(nodes)
    r_i = [nodes[i].r for i in 1:n]

    q_i = [ # Parental instantiations for each node
        prod([r_i[i] for i in inneighbors(graph,node)]) 
        for node in 1:n]

    return [ones(q_i[i], r_i[i]) for i in 1:n]
end

#######################################################
#   Bayesian Score
#######################################################

function bayesian_score_component(M, α)
    p = sum(loggamma.(α + M))
    p -= sum(loggamma.(α))
    p += sum(loggamma.(sum(α,dims=2)))
    p -= sum(loggamma.(sum(α,dims=2) + sum(M,dims=2)))
    return p
end
function bayesian_score(vars, G, D)
    n = length(vars)
    M = statistics(vars, G, D)
    α = prior(vars, G)
    return sum(bayesian_score_component(M[i], α[i]) for i in 1:n)
end

function compute(nodes::Vector{Variable}, samples::Matrix{Int}, starting_graph::DiGraph, restarts = 100)

    best_score = -Inf
    best_graph = nothing
    best_ordering = nothing

    for _ in 1:restarts
        k2_ordering::Vector{Int} = shuffle(1:length(nodes))

        graph = k2_algorithm(k2_ordering, nodes, samples)
        score = bayesian_score(nodes, starting_graph, samples)

        if score > best_score
            best_score = score
            best_graph = graph
            best_ordering = k2_ordering
        end
    end

    # k2_ordering::Vector{Int} = shuffle(1:length(nodes))
        
    # graph = k2_algorithm(k2_ordering, nodes, samples, 20)
    method = LocalDirectedGraphSearch(best_graph, 500)
    graph = local_random_search(method, nodes, samples)
    score = bayesian_score(nodes, graph, samples)
    return graph, score, k2_ordering
end

#######################################################
#   K2 Algorithms
#######################################################

function rand_graph_neighbor(graph, max_nodes::Int)
    n = nv(graph) # Number of vertices
    i = rand(1:n)
    j = mod1(i + rand(2:n)-1, n)
    new_graph = copy(graph)

    # Remove if max edges, or edge already exists
    if outdegree(graph, i) > max_nodes || has_edge(new_graph, i, j)
        rem_edge!(new_graph, i, j)
    else
        add_edge!(new_graph, i, j)

    end

    return new_graph
end

function local_random_search(method::LocalDirectedGraphSearch, nodes::Vector{Variable}, samples::Matrix{Int}, max_nodes::Int = 15)
    graph = method.G
    current_score = bayesian_score(nodes, graph, samples)
    for iter in 1:method.k_max
        # println("Local iteration $iter")
        new_graph = rand_graph_neighbor(graph, max_nodes)
        if is_cyclic(new_graph)
            new_score = -Inf
        else
            new_score = bayesian_score(nodes, new_graph, samples)
        end

        if new_score > current_score
            score, graph = new_score, new_graph
        end
    end

    return graph
end


function k2_algorithm(k2_ordering::Vector{Int}, vars, D, iterations = 100000, max_nodes = 15)
    n = length(vars)

    # Empty graph, to add nodes to
    G = SimpleDiGraph(n)
    for (k,i) in enumerate(k2_ordering[2:end])
        y = bayesian_score(vars, G, D)
        # for iter in 1:iterations
        iter = 0
        while true
            y_best, j_best = -Inf, 0
            for j in k2_ordering[1:k]
            # for x in 1:iterations
                # print("node $x")
                # j = mod1(i + rand(2:n)-1, n)
                
                
                
                if !has_edge(G, j, i)
                    add_edge!(G, j, i)
                    y′ = bayesian_score(vars, G, D)
                    if !is_cyclic(G) && y′ > y_best
                        y_best, j_best = y′, j
                    end
                    rem_edge!(G, j, i)
                end
            end
            
            if y_best > y
                y = y_best
                add_edge!(G, j_best, i)
            else
                break
            end

            iter += 1
            if iter == iterations
                break
            end
        end
        println("Stopped after $iter iterations on node $i: $y")
    end
    return G
end


end # module Project1