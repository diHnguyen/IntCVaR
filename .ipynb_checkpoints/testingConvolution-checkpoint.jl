using JuMP 
using DataFrames
# using Polynomials

include("functionConvolution.jl")

df_cell = DataFrame(CELL = Int64[], Y = Array[], g = Float64[], h = Float64[], LB = Array[], UB = Array[], PROB = Float64[])
push!(df_cell, (1,[1, 1, 1, 1, 1],0.0,0.0,[9, 2, 3, 3, 6],[9, 4, 3, 5, 9],1.0))
push!(df_cell, (2,[0, 1, 1, 0, 1],0.0,0.0,[2, 9, 5, 3, 9],[2, 11, 8, 3, 12],4.0))
push!(df_cell, (3,[1, 1, 1, 1, 0],0.0,0.0,[5, 3, 3, 5, 3],[9, 3, 7, 8, 3],3.0))
# Len = length(cL)
df_cellPoly = DataFrame(CELL = Int64[], df = Any[], NUMPOLY = Int64[], DETSHIFT = Float64[])
# global arcNum = 5
# for k = 1:3
#     global arcNum
#     cL = rand(2:10,arcNum)
#     M = rand(0:4,arcNum)
#     cU = cL +M
#     Y = rand(0:1,arcNum)
#     p = rand(0:5)
#     push!(df_cell, (k,Y,0,0,cL,cU,p))
    
# end
# s = sum(df_cell[k,:PROB] for k = 1:3)
# df_cell[:,:PROB] = df_cell[:,:PROB]./s
println(df_cell)
# global det_Shift = 0
for k = 1:nrow(df_cell)
    cL = df_cell[k,:LB]
    cU = df_cell[k,:UB]
    y = df_cell[k,:Y]
    Len = length(cL)
    println("\nCell ", k)
    unc_arcs = findall((y.>0) .& (cL.<cU))
    det_arcs = findall((y.>0) .& (cL.==cU))
    println("unc_arcs = ", unc_arcs)
    println("det_arcs = ", det_arcs)
    
#     det_arcs = findall((Y.>0) .& (cL.==cU))
#     filter!(x-> !(x in unc_arcs), det_arcs)
#     det_Shift = 0
    if isempty(det_arcs) == false
        det_Shift = sum(cL[i] for i in det_arcs)
    else
        det_Shift = 0
    end
    df, numPoly = Convolve(cL, cU, Len, unc_arcs)
#     df[:,:leftShift] = df[:,:leftShift] .+ det_Shift
    df.w = zeros(numPoly)
    println(df)
    push!(df_cellPoly, (k,df,numPoly, det_Shift))
#     push!(df_cellPoly, (k,DataFrame(),numPoly, det_Shift))
end
# break
# for k = 1:3
# β = 0.7
# α = 20
# gL = 0 #[0,0,0]
# gU = 40 #[20,20,20]
# println(df_cell[:,:PROB])

#FindCVaR(β,α,gL,gU,df_cellPoly, df_cell[:,:PROB])
# println(df_cellPoly)