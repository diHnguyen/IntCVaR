using JuMP 
using Gurobi
using LightGraphs
using DataFrames
using CSV
using TimerOutputs

include("./functionGbound.jl")
include("./function_cornerTest.jl")
include("./myData1.jl")

x_now = zeros(Len)
c_L = copy(cL_orig)
c_U = copy(cU_orig)
c = 0.5*(c_U + c_L)
c_g = c + d.*x_now
y, gx, SP, T, pred, label, path = gx_bound(c_L, c_U, c, c_g, x_now, edge)
println("edge = ", edge)
println("y = ", findall(y.==1))
opt, c_corner = corner_g(y, pred, c_L, c_U, d, x_now, edge, origin, destination)
println(opt)
println(c_corner)