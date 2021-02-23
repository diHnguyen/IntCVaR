using JuMP 
using Gurobi
using LightGraphs
using DataFrames
using CSV
using TimerOutputs

include("./functionGbound.jl")
include("./functionHbound.jl")
include("./functions_SAPath.jl")
include("./myData1.jl")
# println(all_nodes)
# println(cL_orig)

gurobi_env = Gurobi.Env()
epsilon = 1e-3
setparams!(gurobi_env, Heuristics=0.0, Cuts = 0, OutputFlag = 0)

# If we want to add # in Gurobi, then we have to turn of Gurobi's own Cuts 
h1 = Model(() -> Gurobi.Optimizer(gurobi_env))
@variable(h1, 1 >= y_h[1:Len]>=0)
#@variable(h, q[1:Len]>=0)


#Setting constraint for start node
leaving = findall(edge[:,1].== origin)
# println(all_nodes)

@constraint(h1, sum(y_h[k] for k in leaving) == 1)
for i in all_nodes
    if i != destination && i != origin
        incoming = findall(edge[:,2].== i)
        leaving = findall(edge[:,1].== i)
        @constraint(h1, sum(-y_h[k] for k in leaving) + sum(y_h[k] for k in incoming) == 0)
    end
end

# println(Len)
x_now = zeros(Len)
checkedCells = Int64[]
x_count = 0
Delta_Final_Pos = ones(Len).*-1
Delta_Final_Neg = ones(Len).*-1
myLength = length(p)
#                 k=0
#                 while k < myLength #for k = 1:myLength#cell_num in Cell_List
for k = 1:myLength #===SA.1.1===#
#     println("k = ", k)
    if (k in checkedCells) == false #===SA.1.1.1===#
        row = 1 #find(df[:CELL] .== cell_num)[1]
        c_L = copy(cL_orig)
        c_U = copy(cU_orig)
        M = c_U - c_L
        c = 0.5*(c_U + c_L)
        c_g = c + d.*x_now

        r_SA = 0.5*(c_U - c_L)
        
#         c_g = c + d.*x_now

        #FIND TREE, PRED
        y, gx, SP, T, pred, label, path = gx_bound(c_L, c_U, c, c_g, x_now, edge)
        g[k] = gx
        basic, basic_cost, i_basic, nonbasic, nonbasic_cost, i_nonbasic = classify_edges_costs(T,c_g)
        orderedArcList = getSortedBasicArcList(pred, basic, i_basic)

        i_nonbasic_temp = copy(i_nonbasic)
        label_temp = copy(label)
        c_temp = copy(c_g)

        for i = 1:length(i_nonbasic)
            nonbasic_e = i_nonbasic[i]
            u = edge[nonbasic_e,1]
            v = edge[nonbasic_e,2]
            if v in path
                Delta_Final_Neg[nonbasic_e] = label_temp[u] - label_temp[v] + c_temp[nonbasic_e]
            end
        end

        for i = 1 : length(orderedArcList)
            e = orderedArcList[i]
            if y[e] == 1
                myNode = edge[e,2]
                Delta_Final_Pos[e], entering_arc = find_SA_Delta(myNode, Delta_Final_Pos, Delta_Final_Neg, "POS", basic, i_basic, nonbasic, i_nonbasic)
            end
        end
#         println("orderedArcList = ", orderedArcList)
        for i = 1:length(orderedArcList)
            e = orderedArcList[i]
            if Delta_Final_Pos[e] == -1 && Delta_Final_Neg[e] == -1
                myNode = edge[e,2]
                myNonTreeArcList = findall(nonbasic[:,1].== myNode)
                for item in myNonTreeArcList 
                    arc = i_nonbasic[item]
                    ##
                    if Delta_Final_Neg[arc] == -1 #Meaning e has not been checked
                        Delta_Final_Neg[arc] = findNonBasicArc_Loop(label, T, c_g, arc, path, y, pred, basic, i_basic)
                    end
                    ##
                end
                Delta_Final_Neg[e], entering_arc = find_SA_Delta(myNode, Delta_Final_Pos, Delta_Final_Neg, "NEG", basic, i_basic, nonbasic, i_nonbasic)                         
            end
        end

        remainingArcs = findall(Delta_Final_Neg.== -1)
        for e in remainingArcs
            if Delta_Final_Pos[e] == -1 #&& Delta_Final_Neg[e] == -1
                myNode = edge[e,2]
                u = edge[e,1]

                temp_val = label_temp[u] - label_temp[myNode] + c_temp[e]
                Delta_Final_Neg[e], entering_arc = find_SA_Delta(myNode, Delta_Final_Pos, Delta_Final_Neg, "NEG", basic, i_basic, nonbasic, i_nonbasic)                         
                Delta_Final_Neg[e] = Delta_Final_Neg[e]+temp_val
            end
        end

        rule100 = 0
#         println(r_SA)
        rule100 = rule100 + sum(r_SA[e]/Delta_Final_Pos[e] for e = 1:length(Delta_Final_Pos) if Delta_Final_Pos[e] != -1.0) 
        rule100 = rule100 + sum(r_SA[e]/Delta_Final_Neg[e] for e = 1:length(Delta_Final_Neg) if Delta_Final_Neg[e] != -1.0) 
#         println(rule100)
        if rule100 <= 1.0 
            h[k] = g[k] #h[k] = hx
            push!(checkedCells, k)
        end
    end #===End of SA.1.1.1===#
end #===End of SA.1.1===# 
# println(r_SA)
println(Delta_Final_Pos)
println(Delta_Final_Neg)