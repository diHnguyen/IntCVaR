using DataFrames
using JuMP


df_cell = DataFrame(CELL = Int64[], Y = Array[], g = Float64[], h = Float64[], LB = Array[], UB = Array[], PROB = Float64[])
push!(df_cell, (1,[1, 1, 1, 1, 1],0.0,0.0,[9, 2, 3, 3, 6],[9, 4, 3, 5, 9],1.0))
push!(df_cell, (2,[0, 1, 1, 0, 1],0.0,0.0,[2, 9, 5, 3, 9],[2, 11, 8, 3, 12],4.0))
push!(df_cell, (3,[1, 1, 1, 1, 0],0.0,0.0,[5, 3, 3, 5, 3],[9, 3, 7, 8, 3],3.0))

global reps = 10000
global Total = 0
global Count = 0
for k = 1:3
    cL = df_cell[k,:LB]
    cU = df_cell[k,:UB]
    M = cU - cL
    y = df_cell[k,:Y]
#     println(y)
    path = findall(y.==1)
    println(path)
    global Count
    global Total
    global reps
    for r = 1:reps
        pathCost = 0
        for i in path
            myR = rand()*M[i]
            pathCost = pathCost + cL[i] + myR
        end
#         println("path ", pathCost)
        if pathCost <= 21.66015625
            Count = Count + 1 
            Total = Total + pathCost
        end
        
    end
#     println(Total/reps)
    
end

println(Count/reps/3)
println(Total/Count)