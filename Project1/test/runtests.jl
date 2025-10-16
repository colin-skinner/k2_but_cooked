using Project1: Variable, statistics, prior, bayesian_score
using Test
using Graphs

@testset "Statistics Test" begin
    # Sanity check from Example 4.1 of the textbook
    G = SimpleDiGraph(3)
    add_edge!(G, 1, 2)
    add_edge!(G, 3, 2)
    vars = [Variable(:A,2), Variable(:B,2), Variable(:C,2)]
    D = [1 2 2 1; 1 2 2 1; 2 2 2 2]
    M = statistics(vars, G, D)

    @test M[1] == [2 2]
    @test M[2] == [0 0; 0 0; 2 0; 0 2]
    @test M[3] == [0 4]
end



@testset "Prior Test" begin
    # Sanity check from Example 4.1 of the textbook
    G = SimpleDiGraph(3)
    add_edge!(G, 1, 2)
    add_edge!(G, 3, 2)
    vars = [Variable(:A,2), Variable(:B,2), Variable(:C,2)]
    D = [1 2 2 1; 1 2 2 1; 2 2 2 2]
    α = prior(vars, G)

    @test α[1] == [1 1]
    @test α[2] == [1 1; 1 1; 1 1; 1 1]
    @test α[3] == [1 1]
end

@testset "Bayesian Score Test" begin

    nodes = [
        Variable(:parent1, 3)
        Variable(:child1, 3)
        Variable(:parent2, 3)
        Variable(:child2, 3)
        Variable(:parent3, 3)
        Variable(:child3, 3)
    ]

    samples = Matrix([
        3 1 3 1 2 1 1 3 1 3 2 2 2 2 1 2 1 3 3 1
        3 3 1 3 3 3 3 3 3 1 3 3 2 3 1 2 2 1 1 3
        2 2 1 1 3 1 3 2 1 1 3 2 3 3 2 1 1 1 1 1
        3 3 2 2 2 1 1 3 3 3 3 3 1 3 3 3 3 3 3 2
        1 2 3 3 2 1 2 1 2 3 2 3 1 2 3 2 2 1 3 3
        3 3 1 2 1 3 3 3 1 3 3 3 3 3 2 3 3 2 3 1
    ])

    # Graph in old order of nodes
    graph = SimpleDiGraph(length(nodes))
    add_edge!(graph, 1, 2); add_edge!(graph, 1, 4); add_edge!(graph, 3, 4)
    add_edge!(graph, 5, 4); add_edge!(graph, 5, 6)
    
    score = bayesian_score(nodes, graph, samples)

    # https://github.com/sisl/AA228-CS238-Student/blob/main/project1/example/example.score
    @test isapprox(score, -132.57689402451837)
end

# @testset "Bayesian Score Test Different Order" begin

#     nodes = [
#         Variable(:parent1, 3)
#         Variable(:child1, 3)
#         Variable(:parent2, 3)
#         Variable(:child2, 3)
#         Variable(:parent3, 3)
#         Variable(:child3, 3)
#     ]

#     samples = Matrix([
#         3 3 2 3 1 3;
#         1 3 2 3 2 3;
#         3 1 1 2 3 1;
#         1 3 1 2 3 2;
#         2 3 3 2 2 1;
#         1 3 1 1 1 3;
#         1 3 3 1 2 3;
#         3 3 2 3 1 3;
#         1 3 1 3 2 1;
#         3 1 1 3 3 3;
#         2 3 3 3 2 3;
#         2 3 2 3 3 3;
#         2 2 3 1 1 3;
#         2 3 3 3 2 3;
#         1 1 2 3 3 2;
#         2 2 1 3 2 3;
#         1 2 1 3 2 3;
#         3 1 1 3 1 2;
#         3 1 1 3 3 3;
#         1 3 1 2 3 1;
#     ]')
    
#     graph = SimpleDiGraph(length(nodes))
#     add_edge!(graph, 1, 2); add_edge!(graph, 1, 4); add_edge!(graph, 3, 4)
#     add_edge!(graph, 5, 4); add_edge!(graph, 5, 6)

#     order = [1 4 2 5 3 6]

#     ordered_samples = similar(samples, Int)

#     for i in order
#         ordered_samples[i,:] = samples[order[i],:]
#     end

#     score = bayesian_score(nodes, graph, ordered_samples)

#     # https://github.com/sisl/AA228-CS238-Student/blob/main/project1/example/example.score
#     @test isapprox(score, -132.57689402451837)
# end