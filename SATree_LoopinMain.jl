using JuMP 
using Gurobi
using LightGraphs
using DataFrames
using CSV
using TimerOutputs

include("./functionGbound.jl")
include("./functionHbound.jl")
include("./functions_SATree.jl")
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
    if (k in checkedCells) == false #===SA.1.1.1===#
        row = 1 #find(df[:CELL] .== cell_num)[1]
        c_L = copy(cL_orig)
        c_U = copy(cU_orig)
        M = c_U - c_L
        c = 0.5*(c_U + c_L)
        c_g = c + d.*x_now

        r_SA = 0.5*(c_U - c_L)

        y, gx, SP, T, pred, label, path = gx_bound(c_L, c_U, c, c_g, x_now, edge)

        g[k] = gx
        basic, basic_cost, i_basic, nonbasic, nonbasic_cost, i_nonbasic = classify_edges_costs(T,c_g)
        orderedArcList = getSortedBasicArcList(pred, basic, i_basic)

        label_temp = copy(label)
        c_temp = copy(c_g)

        for e in i_nonbasic
            u = edge[e,1]
            v = edge[e,2]
            Delta_Final_Neg[e] = label_temp[u] - label_temp[v] + c_temp[e]
        end

        for i = 1 : length(orderedArcList)
            e = orderedArcList[i]
            myNode = edge[e,2]
            if y[e] == 1 #if arc is on the shortest path y
                Delta_Final_Pos[e] = find_basic_Delta(myNode, Delta_Final_Pos, Delta_Final_Neg, "POS", basic, i_basic, nonbasic, i_nonbasic)
                #Do not have to find Delta_Final_Neg : all set to (-1)
            else
                Delta_Final_Neg[e] = find_basic_Delta(myNode, Delta_Final_Pos, Delta_Final_Neg, "NEG", basic, i_basic, nonbasic, i_nonbasic)
            end
        end
        rule100 = 0
        rule100 = rule100 + sum(r_SA[e]/Delta_Final_Pos[e] for e = 1:length(Delta_Final_Pos) if Delta_Final_Pos[e] != -1.0) 
        rule100 = rule100 + sum(r_SA[e]/Delta_Final_Neg[e] for e = 1:length(Delta_Final_Neg) if Delta_Final_Neg[e] != -1.0) 
        if rule100 <= 1.0 
            h[k] = g[k] #h[k] = hx
            push!(checkedCells, k)
        end

    end #if (k in checkedCells) == false
        #################################SENSITIVITY ENDS#################################
end 

println(Delta_Final_Pos)
println(Delta_Final_Neg)

# y, hx = hx_bound(h, c_L, c_U, c, d, x_now, edge)
# Delta_Final = zeros(Len)
# cell_num = 0
# Cell_List = [1]
# myLength = length(Cell_List)
# global k=0

# while k < myLength #for k = 1:myLength#cell_num in Cell_List
#     k=k+1
#     cell_num = Cell_List[k]

#     #k = Cell_List[cell_num]
#     row = 1 #find(df[:CELL] .== cell_num)[1]
#     c_L = copy(cL_orig)
#     c_U = copy(cU_orig)
#     c = 0.5*(c_U + c_L)
#     c_g = c + d.*x_now
    
#     r_SA = 0.5*(c_U - c_L)
#     c = (c_U + c_L)/2
#     M = c_U - c_L
#     c_g = c + d.*x_now

#     #FIND TREE, PRED
#     y, gx, SP, T, pred, label, path = gx_bound(c_L, c_U, c, c_g, x_now, edge)
#     g[cell_num] = gx
#     println("1")
#     ###println(f,"y = ", y)
#     ###println(f,"T = ", T)
# #                 println(f,"c_U = ", c_U)
# #                 println(f,"c_L = ", c_L)  
#     for e = 1:Len

        
#         if r_SA[e] > 0
#             println("Acr #", e)
#             #taboo_list = []
#             basic, basic_cost, i_basic, nonbasic, nonbasic_cost, i_nonbasic = classify_edges_costs(T,c_g)
#             ##println("c = ", c_g)
#             label_temp = copy(label)
#             Δ = 0
#             δ = []
#             δ_ind = 0
#             T_temp = copy(T)
#             pred_temp = copy(pred)
#             c_temp = copy(c_g)
#             edge_num = copy(e)
#             N = Int[]
#             cur_edge = Int[]
#             path_temp = copy(path)
#             arc_type = what_arc(T_temp, y, e)
#             new_SP = false
#             println("Type: ", arc_type)
#             if arc_type == 2
#                 u = edge[e,1]
#                 v = edge[e,2]
#                 t = pred[v]
#                 e_temp = [t, v]
# #                 println("A")
#                 ###println(f,"(t v) = (",t, " ", v, ")")
# #                 B = findall((basic .== e_temp,2)) #find a tree edge to be replaced by nonbasic e
# #                 i = B[1]
#                 ###println(f, "edge to be replaced = ", i)
# #                 i = i_basic[i]
#                 i = 0
#                 println(e_temp)
#                 for temp_i in i_basic
#                     println(edge[temp_i,:])
#                     if edge[temp_i,:] == e_temp
#                         i = temp_i
#                     end
#                 end
# #                 println("B")

#                 val = label_temp[u] - label_temp[v] + c_temp[e]
#                 push!(δ, val)
#                 ###println(f,"δ = ", δ)
#                 δ_min = minimum(δ)
#                 c_temp[e] = c_temp[e] - δ_min
#                 Δ = Δ + δ_min
#                 pred_temp[v] = u
#                 ###println(f, "Tree before replace: ", T_temp)
#                 T_temp[e] = 1
#                 T_temp[i] = 0

#                 ###println(f, "Tree after replace: ", T_temp)
#                 if T_temp == T
#                     new_SP = false  
#                 else
#                     new_SP = true
#                     Delta_Final[e] = Δ
#                 end
#                 arc_type = what_arc(T_temp, y, e)
#                 ###println(f, "new_SP = ", new_SP)
#             end
            
#             println("Before new_SP Loop")
#             loop = 0
#             println("new_SP = ", new_SP)
#             while new_SP == false
# #                 println()
#                 loop = loop+1
#                 println("Loop ", loop)
#                 basic, basic_cost, i_basic, nonbasic, nonbasic_cost, i_nonbasic = classify_edges_costs(T_temp,c_temp)
                
#                 ##println( "Current T = ", T_temp)
#                 ##println( "Current y = ", y)
#                 ##println("Current label = ", label_temp)

#                 if arc_type == 0 #THIS STATEMENT IS WORKING
#                     N = find_affected_node(T_temp, pred_temp, e, edge)
#                     #Hereeeeee: Must update tree, path
#                     ##println("N = ", N)
#                     δ_min, δ_ind, temp = SA_basicSP(N, nonbasic, nonbasic_cost, label_temp, c_temp[e])
#                     ##println("Δ = ", Δ)
#                     ##println("δ_min = ", δ_min)
#                     if δ_min < 0
#                         Delta_Final[e] = δ_min
#                         new_SP = true
#                     else
#                         ##println("nonbasic arcs: ", nonbasic)
#                         ##println("index of arc moving into the basis: ", δ_ind)
#                         u = nonbasic[δ_ind, 1]
#                         v = nonbasic[δ_ind, 2]
#                         t = pred_temp[v]
#                         ##println("pred = ", pred_temp)
#                         ##println("u = ", u)
#                         ##println("v = ", v)
#                         ##println("t = ", t)
#                         e_temp = [t v]
#                         ##println("Basic arc becoming nonbasic: ", e_temp)
#                         #push!(taboo_list, i_nonbasic[δ_ind])
#                         ##println("Nonbasic arc becoming basic: ", [u v])
#                         ##println("basic arcs = ", basic)
#                         B = find(all(basic .== e_temp,2)) #find a tree edge to be replaced by nonbasic e
#                         ##println("B = ", B)
#                         i = B[1]
#                         i = i_basic[i]
#                         pred_temp[v] = nonbasic[δ_ind, 1]
#                         c_temp[e] = temp
#                         T_temp[i_nonbasic[δ_ind]] = 1
#                         T_temp[i] = 0 
#                         Δ = Δ + δ_min
#                         Delta_Final[e] = Δ
#                         ##println("CHECK MY CURRENT TREE: ", T_temp)
#                     end
#                 end
#                 if arc_type == 1
#                     N = find_affected_node(T_temp, pred_temp, e, edge)
#                     ##println("N = ", N)
#                     δ_min, δ_ind, temp = SA_basic(N, nonbasic, nonbasic_cost, label_temp, c_temp[e])
#                     ##println("Δ = ", Δ)
#                     ##println("δ_min = ", δ_min)
#                     if δ_min < 0
#                         Delta_Final[e] = δ_min
#                         new_SP = true
#                     else
#                         ##println("nonbasic arcs: ", nonbasic)
#                         ##println("index of arc moving into the basis: ", δ_ind)
#                         u = nonbasic[δ_ind, 1]
#                         v = nonbasic[δ_ind, 2]
#                         t = pred_temp[v]
#                         e_temp = [t v]
#                         ##println("Basic arc becoming nonbasic: ", e_temp)
#                         ##println("Nonbasic arc becoming basic: ", [u v])
#                         B = find(all(basic .== e_temp,2)) #find a tree edge to be replaced by nonbasic e
#                         i = B[1]
#                         i = i_basic[i]
#                         pred_temp[v] = nonbasic[δ_ind, 1]
#                         c_temp[e] = temp
#                         T_temp[i_nonbasic[δ_ind]] = 1
#                         T_temp[i] = 0 ################
#                         Δ = Δ + δ_min
#                         Delta_Final[e] = Δ
#                         ##println("Pred = ", pred_temp)
#                     end

#                 end 
#                 #CHECK FOR LOOP:

#                 println("Before T_temp IF")
#                 if T_temp != T
#                     new_SP = true
#                 else
#                     loop = false
#                     if pred_temp[origin] != 0
#                         loop = true
#                     else
#                         v = edge[e,2]
#                         visit = collect(1:last_node)
#                         ##println("visit = ", visit)
#                         while v != origin && visit[v] != 0
#                             ##println("v = ", v)
#                             ##println("pred[v] = ", pred_temp[v] )
#                             visit[v] = 0
#                             v = pred_temp[v]

#                             ##println("visit = ", visit)

#                         end

#                         if v != origin
#                             loop = true
#                             ##println("HERE")
#                         end
#                     end
#                     if loop == true
#                         new_SP = true

#                     end
#                 end
#                 ##println("new_SP = ", new_SP)

#             end
#         else
#             Delta_Final[e] = -1.0
#         end
#     end
#     ###println(f,"Deltas = ", Delta_Final)
#     rule100 = 0

#     rule100 = sum(r_SA[e]/Delta_Final[e] for e = 1:length(Delta_Final) if Delta_Final[e] != -1.0) 

#     #println(f,"cell ", cell_num)
# #                 println(f,"100% Rule Result for cell ", cell_num, ": ", rule100)

#     ###println(f,"r_SA = ", r_SA)
#     ###println(f,"Delta_Final = ", Delta_Final)
#     println("rule100")
#     if rule100 <= 1.0 
#         println("cell_num = ", cell_num)
#         y_h, hx = hx_bound(c_L, c_U, d, x_now) 
#         println("A")
#         h[cell_num] = hx #h[k] = hx
#         Cell_List[k], Cell_List[myLength] = Cell_List[myLength], Cell_List[k]
#         pop!(Cell_List)
#         println("B")
#         k = k - 1
#         myLength = myLength - 1

#         y_h, hx = hx_bound(c_L, c_U, d, x_now)
        
#         if length(h) >= k #cell_num
#             h[cell_num] = hx 
#             #h[k] = hx
#         else
#             push!(h,hx)
#             print(1)
#         end
# #                     println(f, "g = ", g)
# #                     println(f, "h = ", h)
#     end

#     ######REMOVE THIS AFTER DONE TESTING FOR SA    
# #                 println(f, "g = ", g[cell_num])
#     y_h, hx = hx_bound(c_L, c_U, d, x_now)
#     h[cell_num] = hx 
#     stopping_cond = delta2 + 1
# end

# println(Delta_Final)




