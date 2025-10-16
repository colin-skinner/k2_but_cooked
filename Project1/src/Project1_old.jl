module Project1

using Graphs
using Printf
using DataStructures

"""
    import_data(filename::String, debug::Bool = true)::Tuple{Array{String}, Array{Array{Int8}}}

Takes in a filename and outputs an array of node names and a nested array of samples and their values.

e.g.
```
SubString{String}["parent1", "child1", "parent2", "child2", "parent3", "child3"]
Int8[3, 3, 2, 3, 1, 3]
Int8[1, 3, 2, 3, 2, 3]
...
Int8[1, 3, 1, 2, 3, 1]
```
"""
function import_data(filename::String, debug::Bool = false)::Tuple{Array{String}, Vector{Vector{Int8}}}

    node_names = String[]
    samples = Vector{Vector{Int8}}()

    open(filename, "r") do io

        # Break if the file is empty
        header = readline(io)

        if length(header) > 0

            # Parse node names
            node_names = split(header, ',')
            map!(String, node_names)

            if debug
                @printf("%s\n", node_names)
            end

            # Parse samples
            for line in eachline(io)
                values_strings = split(line, ',')
                data = parse.(Int8, values_strings)
                push!(samples, data)
                
                if debug
                    @printf("%s\n", data)
                end

            end
        end
        
    end

    
    return node_names, samples
end

"""
    write_gph(dag::DiGraph, idx2names, filename)

Takes a DiGraph, a Dict of index to names and a output filename to write the graph in `gph` format.
"""
function write_gph(dag::DiGraph, idx2names::Array{String}, filename)
    open(filename, "w") do io
        for edge in edges(dag)
            @printf(io, "%s,%s\n", idx2names[src(edge)], idx2names[dst(edge)])
        end
    end
end


function compute(node_names::Array{String}, samples::Vector{Vector{Int8}}, outfile::String)

    n = length(node_names)


    node_num_to_name = [
            "parent1",
            "parent2",
            "parent3",
            "child1",
            "child2",
            "child3"
        ]

    # Creating graph (to be done each iteration)
    graph = SimpleDiGraph(n)
    add_edge!(graph, 1, 5); add_edge!(graph, 2, 5); add_edge!(graph, 3, 5)
    add_edge!(graph, 1, 4); add_edge!(graph, 3, 6)
    
    # Calculating the score

    order = [1,2,3,4,5,6]
    score = bayesian_score(graph, samples, order, n)
    @printf("Score: %d\n", score)
    
    # Writing results
    write_gph(graph, node_num_to_name, outfile)
    # WRITE YOUR CODE HERE
    # FEEL FREE TO CHANGE ANYTHING ANYWHERE IN THE CODE
    # THIS INCLUDES CHANGING THE FUNCTION NAMES, MAKING THE CODE MODULAR, BASICALLY ANYTHING

end

function bayesian_score(graph::SimpleDiGraph, samples::Vector{Vector{Int8}}, order::Array{Int64}, n::Int64)

    # Allocating space for all of the m_indices (undefined values so far)
    m_array = Array{Array{Int8}}(undef, n)

    # HORRIBLY SPACE INEFFICIENT. FIX WITH MAP?
    for i in 1:n
        m_array[i] = zeros(Int8, n)
    end

    for sample in samples
        # Looping through in the topological order
        for i in order
            parents = inneighbors(graph, i)
            # print("$parents\n")
        
            
            num_parents = length(parents)
            # parent_instantiations = 3^length(num_parents)
            
            # print("For x_$i, out of $parent_instantiations\n")

            counter = num_parents
            
            # ASSUMES 3 EACH FOR THE EXAMPLE. CHANGE LATER
            NUM_PARENTAL_VALUES = 3
            j = 1
            # Look up parent values
            for parent in parents
                # DID MATH UNTIL IT WORKED. CHANGE ORDER TO MAKE IT BETTER
                j += (NUM_PARENTAL_VALUES^(counter - 1) * (sample[parent] - 1))

                counter -= 1
            end

            k = sample[i]
            print("$i $j $k\n") # GOOD FOR DEBUGGING

            
            print("\n")

        end

        break
            

    end

    
    return 0
end

if !(length(ARGS) in [2,3])
    error("usage:\n\tjulia project1.jl <infile>.csv <outfile>.gph\n\tjulia project1.jl <infile>.csv <outfile>.gph debug[bool]")
end

inputfilename::String = ARGS[1]
outputfilename::String = ARGS[2]

# Maybe put -d as a flag
if length(ARGS) == 3
    debug::Bool = parse(Bool, ARGS[3])
    node_names, samples = import_data(inputfilename, debug)
else
    node_names, samples = import_data(inputfilename)
end




compute(node_names, samples, outputfilename)


end