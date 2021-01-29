#aqui vai estar o modelo!!!!

import Pkg
using JuMP
using Gurobi

# data struct to store the layout
mutable struct CutDataQ
    t #number of cut edges
    weights #weight of each edge
    C::Array{Pair} #list of cut edges
end

function readDataQ(f)
    d = [[]]
    w = [[]]
    x = []
    y = []
    C = []
    v = 0
    n = 0
    i = 1
    t = 0

    #reading each file of the file
    for l in eachline(f)

        #splitting the line by "space"
        q = split(l, " ")

        a = q[1]
        if a == "p" #number of vertices
            n = parse(Int64, q[2])
            d = zeros(n,n)
            x = [0.0 for i = 1:n]
            y = [0.0 for i = 1:n]
            w = [0.0 for i = 1:n, j = 1:n]
        elseif a == "v" # reads a vertex
            v = parse(Int64,q[2])
            x[v] = parse(Float64,q[3])
            y[v] = parse(Float64,q[4])
        elseif a == "b" #number of edges
            sz = parse(Int64,q[2])
            w = zeros(2*sz,2*sz)
            C = [Pair(0,0) for i = 1:sz]
            t = sz
            i = 1
        elseif a == "e" #reads an edge
            C[i] = Pair(parse(Int64,q[2]),parse(Int64,q[3]))
            i = i+1
        else
        end
    end

    maxd = 0
    #compute the distances between the vertices
    for i = 1:n
        for j = 1:n
            d[i,j] = sqrt((abs(x[i]- x[j]))^2 + (abs(y[i] - y[j]))^2)
            if maxd < d[i,j]
                maxd = d[i,j]
            end
        end
    end
    maxd = maxd + 1

    #edge - edge
    for i = 1:t
        for j = 1:t
            v1 = C[i].first
            u1 = C[i].second
            v2 = C[j].first
            u2 = C[j].second

            if u1 != v2
                w[i,j] = d[u1,v2]
            else
                w[i,j] = 0
            end

            if u1 != u2
                w[i,j+t] = d[u1,u2]
            else
                w[i,j+t] = 0
            end

            if v1 != v2
                w[i+t,j] = d[v1,v2]
            else
                w[i+t,j] = 0
            end

            if v1 != u2
                w[i+t,j+t] = d[v1,u2]
            else
                w[i+t,j+t] = 0
            end

        end
    end

    for i in 1:t
        w[i,i] = maxd
        w[i+t,i] = maxd
        w[i, i+t] = maxd
        w[i+t,i+t] = maxd
    end

    return CutDataQ(t,w,C)
end


function QuadraticFormulation(data::CutDataQ)
    #Using CPLEX
    md = Model(solver = GurobiSolver())
    t = data.t
    w = data.weights
    C = data.C

    #creating the variables
    @variable(md, x[i=1:2*t , j=1:t ] , Bin)

    for i in 1:t
        @constraint(md, sum(x[e, i] for e = 1:2*t) == 1)
    end

    for e in 1:t
        @constraint(md, sum(x[e, i] + x[e+t,i] for i = 1:t) == 1)
    end

    @objective(md, Min, sum(w[i,k]*x[i,j]*x[k,j+1] + w[i+t,k]*x[i+t,j]*x[k,j+1] + w[i,k+t]*x[i,j]*x[k+t,j+1] for i = 1:t, k = 1:t , j=1:t-1 ) )

    return md
end


function LinearizationFormulation(data::CutDataQ)
    #Using CPLEX
    md = Model(solver = GurobiSolver())
    t = data.t
    w = data.weights
    C = data.C

    #creating the variables
    @variable(md, x[i=1:2*t , j=1:t ] , Bin)
    @variable(md, y[i=1:2*t , j=1:2*t , k=1:t ] , Bin)

    for i in 1:t
        @constraint(md, sum(x[e, i] for e = 1:2*t) == 1)
    end

    for e in 1:t
        @constraint(md, sum(x[e, i] + x[e+t,i] for i = 1:t) == 1)
    end

    for e in 1:2*t
        for l in 1:2*t
            for i in 1:(t-1)
                @constraint(md, y[e,l,i] >= x[e,i] + x[l,i+1] - 1)
            end
        end
    end

    @objective(md, Min, sum(w[i,k]*y[i,k,j] + w[i+t,k]*y[i+t,k,j] + w[i,k+t]*y[i,k+t,j] + w[i+t,k+t]*y[i+t,k+t,j] for i = 1:t, k = 1:t , j=1:t-1 ) )

    return md
end


function QuadraticHeuristic(data::CutDataQ, α)
    #QUADRATIC MODEL
    md = Model(solver = GurobiSolver())
    t = data.t
    w = data.weights
    C = data.C
    #creating the variables
    @variable(md, 0 <= x[i=1:2*t , j=1:t ] <= 1)
    for i in 1:t
        @constraint(md, sum(x[e, i] for e = 1:2*t) == 1)
    end
    for e in 1:t
        @constraint(md, sum(x[e, i] + x[e+t,i] for i = 1:t) == 1)
    end
    @objective(md, Min, sum(w[i,k]*x[i,j]*x[k,j+1] + w[i+t,k]*x[i+t,j]*x[k,j+1] + w[i,k+t]*x[i,j]*x[k+t,j+1] for i = 1:t, k = 1:t , j=1:t-1 ) )
    solve(md)

    x_val = values(x)

    Δ = [2*t for i in 1:2*t]
    Γ = Set([i for i in 1:2*t])
    L = Set([i for i in 1:t])
    T = [0 for i in 1:t]
    β = [-1 for i in 1:2*t]
    size = α*t

    for id in 1:size
        for i in Γ
            T[i] = 1
            for j in L
                if x_val[i, j] >= x_val[i, T[i]]
                    x_val[i, T[i]] = x_val[i, j]
                    T[i] = j
                end
            end
        end


        for i in Γ
            for j in L
                if j != T[i] && Δ[i] > x_val[i, T[i]] - x_val[i, j]
                    Δ[i] = x_val[i, T[i]] - x_val[i, j]
                end
            end
        end

        Ð = 0
        γ = 1
        for i in Γ
            if Ð < Δ[i]
                Ð = Δ[i]
                γ = i
            end
        end


        β[γ] = T[γ]
        setdiff(Γ, [γ])
        setdiff(L, [T[γ]])
    end

    md_int = Model(solver = GurobiSolver())
    #creating the variables
    @variable(md_int, x_int[i=1:2*t , j=1:t ], Bin)
    for i in 1:t
        @constraint(md_int, sum(x_int[e, i] for e = 1:2*t) == 1)
    end

    for e in 1:t
        @constraint(md_int, sum(x_int[e, i] + x_int[e+t,i] for i = 1:t) == 1)
    end

    for i in 1:2*t
        if β[i] > 0
            @constraint(md_int, x_int[i, B[i]] == 1)
        end
    end

    @objective(md_int, Min, sum(w[i,k]*x_int[i,j]*x_int[k,j+1] + w[i+t,k]*x_int[i+t,j]*x_int[k,j+1] + w[i,k+t]*x_int[i,j]*x_int[k+t,j+1] for i = 1:t, k = 1:t , j=1:t-1 ) )

    return md_int

end



f = open("/home/marcio/Dropbox/WORK/code/Cortes-Julia/grafo.col")
h = readDataQ(f)

println("OK!")

#mdl = LinearizationFormulation(h)

#print(mdl)

#status = solve(mdl)
#println("status = ", status)

#bestsol = getobjectivevalue(mdl)
#println("bestsol = ", bestsol)

#bestbound = getobjbound(mdl)
#println("bestbound = ", bestbound)

#solvernodes = getnodecount(mdl)
#println("Nodes = ", solvernodes)

#solvertime = getsolvetime(mdl)
#println("Solver time = ",solvertime)



#mdll = QuadraticFormulation(h)

#print(mdll)

#status = solve(mdll)
#println("status = ", status)

#bestsol = getobjectivevalue(mdll)
#println("bestsol = ", bestsol)

#bestbound = getobjbound(mdll)
#println("bestbound = ", bestbound)

#solvernodes = getnodecount(mdll)
#println("Nodes = ", solvernodes)

#solvertime = getsolvetime(mdll)
#println("Solver time = ",solvertime)


mdh = QuadraticHeuristic(h, 0.5)
