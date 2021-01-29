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



function updateLambda(n, λ, x, C::Array{Pair}, π)

    for i = 1:length(C)
        λ[i] = λ[i] + π*(1 - 1*(x[C[i].first , C[i].second]) - 1*(x[C[i].second, C[i].first]))
    end
end

function updateWeights(n, λ, w, weights,C::Array{Pair})
  for i = 1:n
    for j = 1:n
      w[i,j] = weights[i,j]
    end
  end
  for i = 1:length(C)
      w[C[i].first, C[i].second] = w[C[i].first, C[i].second] - λ[i]
      w[C[i].second, C[i].first] = w[C[i].first, C[i].second] - λ[i]
  end
end


function lagrangian(data::CutData)
  wg = data.weights
  n = data.n
  C = data.C
  w = [0.0 for i = 1:n, j=1:n]
  for i =1:n
    for j = 1:n
      w[i,j] = wg[i,j]
    end
  end
  λ = [0 for i=1:n]
  best_λ = λ
  π = 2
  ϵ = 0.01
  sol = -typemax(Float64)
  sol_c = 0
  imp = 0
  p = nothing

  md = Model(solver = GurobiSolver())
  @variable(md, x[i=1:n , j=1:n ] , Bin)
  @variable(md, y[i=1:n] , Bin)

  @constraint(md, sum(x[1, v] for v = 1:n) == 1)

  @constraint(md, sum(y[v] for v = 1:n) == 1)

  for v in 2:n
      @constraint(md, sum(x[u, v] for u = 1:n) - sum(x[v, u] for u = 1:n) - y[v] == 0)
  end

  val_x = [0.0 for i = 1:n, j=1:n]
  it = 1

  while (imp < 3) & (it < 3)
    it = it+1

    println("UPDATING WEIGHTS")
    updateWeights(n,λ, w,wg, C)
    println("WEIGHTS UPDATED WITH SUCESS")

    @objective(md, Min, sum(w[i,j]*x[i,j] for i = 1:n, j=1:n))

    println("SOLVING MODEL")

    print(md)

    solve(md)

    println("OPTIMIZATION OK")

    sol_c = getobjectivevalue(md)

    println("OPTIMAL SOLUTION:")
    println(sol_c)

    for i =1:n
      for j =1:n
        val_x[i,j] = getvalue(x[i,j])
      end
    end

    println("UPDATING LAMBDA")
    updateLambda(n,λ, val_x, C, π)
    println("UPDATE LAMBDA OK")


    if (sol_c - sol) > ϵ
      imp = imp + 1
      sol = sol_c
      println("IMPROVEMENT FOUND, UPDATING DATA")
    else
      imp = 0
      println("NO IMPROVEMENT FOUND")
    end

  end

  lb = sol

  return lb
end

f = open("/home/marcio/Dropbox/WORK/code/Cortes-Julia/grafo.col")
h = readData(f)

print(h)

println("OK!")

sl = lagrangian(h)

println("BEST SOLUTION:")

print(sl)
