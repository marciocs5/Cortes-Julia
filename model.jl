#aqui vai estar o modelo!!!!

import Pkg
using JuMP
using Gurobi

# data struct to store the layout
mutable struct CutData
    n #number of vertices
    weights #weight of each edge
    C::Array{Pair} #list of cut edges
end

function readData(f)
    w = [[]]
    x = []
    y = []
    C = []
    v = 0
    n = 0
    i = 1

    #reading each file of the file
    for l in eachline(f)

        #splitting the line by "space"
        q = split(l, " ")

        a = q[1]
        if a == "p" #number of vertices
            n = parse(Int64, q[2])
            x = [0.0 for i = 1:n]
            y = [0.0 for i = 1:n]
            w = [0.0 for i = 1:n, j = 1:n]
        elseif a == "v" # reads a vertex
            v = parse(Int64,q[2])
            x[v] = parse(Float64,q[3])
            y[v] = parse(Float64,q[4])
        elseif a == "b" #number of edges
            sz = parse(Int64,q[2])
            C = [Pair(0,0) for i = 1:sz]
            i = 1
        elseif a == "e" #reads an edge
            C[i] = Pair(parse(Int64,q[2]),parse(Int64,q[3]))
            i = i+1
        else
        end
    end

    #compute the distances between the vertices
    for i = 1:n
        for j = 1:n
            w[i,j] = sqrt((abs(x[i]- x[j]))^2 + (abs(y[i] - y[j]))^2)
        end
    end
    return CutData(n,w,C)
end



function FlowFormulation(data::CutData)
    #Using CPLEX
    md = Model(solver = GurobiSolver())
    n = data.n
    w = data.weights
    C = data.C

    #creating the variables
    @variable(md, x[i=1:n , j=1:n ] , Bin)
    @variable(md, y[i=1:n] , Bin)

    #sum of all elements
    @constraint(md, sum(x[1, v] for v = 1:n) == 1)

    @constraint(md, sum(y[v] for v = 1:n) == 1)

    for v in 2:n
        @constraint(md, sum(x[u, v] for u = 1:n) - sum(x[v, u] for u = 1:n) - y[v] == 0)
    end

    for (u,v) in C
        @constraint(md, x[u, v] + x[v,u] == 1)
    end

    @objective(md, Min, sum(w[i,j]*x[i,j] for i = 1:n, j=1:n))

    return md
end


f = open("/home/marcio/Dropbox/WORK/code/Cortes-Julia/grafo.col")
h = readData(f)

print(h)

println("OK!")

mdl = FlowFormulation(h)

print(mdl)


status = solve(mdl)
println("status = ", status)

bestsol = getobjectivevalue(mdl)
println("bestsol = ", bestsol)

bestbound = getobjbound(mdl)
println("bestbound = ", bestbound)

#opengap = getobjgap(mdl)
#println("Nodes = ", opengap)

solvernodes = getnodecount(mdl)
println("Nodes = ", solvernodes)

solvertime = getsolvetime(mdl)
println("Solver time = ",solvertime)
