#Ready for upload
using JuMP 
using Gurobi
using LightGraphs
using DataFrames, Query
using CSV
using TimerOutputs
using Dates
using Polynomials

myRun = Dates.format(now(), "HH:MM:SS")
global gurobi_env = Gurobi.Env()
global edge, cL_orig, cU_orig, Len, c_orig, yy, SP_init, p,g,h, origin, destination, last_node, all_nodes, M_orig, delta1, delta2, b, last_node 
global β = 0.01
# gurobi_env.setParam("LogToConsole", 0)

to = TimerOutput()
# myFile = "./myData1.jl"
myFile = "./Instances_Paper1/Center Instances/20NodesCenter_1.jl"
include(myFile)
include("PartitionRules.jl")
include("functionGbound.jl")
include("functionHbound.jl")
include("functionPartition_BasicAlg.jl")
include("functionConvolution.jl")
include("functionFindCVaR.jl")

println("Running...")
global epsilon = 1e-4
global delta3 = 1.0
setparams!(gurobi_env, Heuristics=0.0, Cuts = 0, OutputFlag = 0)

# If we want to add # in Gurobi, then we have to turn of Gurobi's own Cuts 
h1 = Model(() -> Gurobi.Optimizer(gurobi_env))
# h1.setParam("OutputFlag", 0)
# set_optimizer_attribute(h1, "OutputFlag", 0)

@variable(h1, 1 >= y_h[1:Len]>=0)
#@variable(h, q[1:Len]>=0)


#Setting constraint for start node
leaving = findall(edge[:,1].== origin)


#Setting constraints for remaining none-sink/start nodes
@constraint(h1, sum(y_h[k] for k in leaving) == 1)
for i in all_nodes
    if i != destination && i != origin
        incoming = findall(edge[:,2].== i)
        leaving = findall(edge[:,1].== i)
        @constraint(h1, sum(-y_h[k] for k in leaving) + sum(y_h[k] for k in incoming) == 0)
    end
end

#MAIN PROGRAM:
df_constraints = DataFrame(NUM = Int[], CELL = Int[], Y = Array[], SP = Float64[])
df_cell = DataFrame(CELL = Int[], Y = Array[], Y_Lk = Array[], g = Float64[], h = Float64[], gL = Float64[], LB = Array[], UB = Array[], PROB = Float64[])


MP_obj = 0.0

push!(df_cell, (1, yy,yy, SP_init, 0, 0, cL_orig, cU_orig, 1))
push!(df_constraints, (1, 1,yy,SP_init))
##println(f,"MASTER PROBLEM==========================================================================================")

zNum = 5
cRefNum = 20
m = Model(() -> Gurobi.Optimizer(gurobi_env)) # If we want to add # in Gurobi, then we have to turn of 
# set_optimizer_attribute(m, "OutputFlag", 0)    #Gurobi's own Cuts 
# println("1")
@variable(m, x[1:Len], Bin)
@variable(m, α)
@variable(m, 1e6 >= z[1:zNum] >= 0)
# @constraintref constr[1:200000]
# @ConstrRef constr[1:200000]

# constr = Array{JuMP.JuMPArray{JuMP.ConstraintRef,1,Tuple{Array{Int64,1}}}}()

constr = Array{JuMP.ConstraintRef}(undef, cRefNum)
@constraint(m, sum(x[i] for i=1:Len) == b) #attack budget = 2
constr[1] = @constraint(m, α <= SP_init + sum(yy[i]*x[i]*d[i] for i=1:Len) + z[1]  )
@objective(m, Max, α - (1/(1-β))*sum(p[i]*z[i] for i = 1:length(p)) )#w - sum(s[k] for k=1:length(s))/length(s) )
global x_sol = []
global α_sol = 0
global z_sol = []
global x_now = []
global α_now = 0
global z_now = []
global last_x = zeros(Len)
global con_num = 1
global newCell = 1
global total_time = 0.0
global iter = 0
# global K = Int64[1]
global K_bar = Int64[1]
global LB = 0
global MP_obj = 1e6
global K_newly_added = []
global K_removed = []
start = time()
global terminate_cond = false

while terminate_cond == false #&& iter < 5#&& isempty(K_bar) == false
    global α, β, iter, total_time, K_bar, K_newly_added, K_removed, LB, MP_obj, con_num, newCell #, min_gLk, max_gUk
    global x_sol, z_sol, α_sol, last_x, x_now,α_now,z_now, terminate_cond
#     println("lengthK = ", length(K))
    while isempty(K_bar) == false #length(K_bar) > K
#         stopping_cond = 0
        iter = iter + 1
        
        optimize!(m) 
        println("\nIter : ", iter," ; LB = ", LB)
#         println(m)
#         println("", df_cell)
#         K = vcat(K, K_newly_added)
        
        if termination_status(m) == MOI.OPTIMAL
            MP_obj = JuMP.objective_value.(m)
            x_now = JuMP.value.(x)
            α_now = JuMP.value.(α)
#             last_x = x_now
            z_now = JuMP.value.(z)
            println("Obj = ", MP_obj, "\tα_now = ", α_now)
            println("z_now = ", z_now[1:nrow(df_cell)])
            println("Interdiction ", findall(x_now.==1))
            for i in eachrow(df_cell)
               println("Cell ", i.CELL, " : ", i.g, "\t", findall(i.Y .>0)) 
            end
            if iter == 5
               println(m) 
            end
        end
        
        if termination_status(m) != MOI.OPTIMAL || MP_obj <= LB
            terminate_cond = true
            K_bar = []
            if termination_status(m) != MOI.OPTIMAL
                println("MP not feasible...")
            else
                println("MP_Obj did not improve..." )
            end
        else
            O1Flag = true
            if last_x != x_now
                K_bar = collect(1:newCell)
                last_x = x_now
                for k in K_bar #_partition  
#                     row = findall(df_cell[!,:CELL] .== k)[1]
                    c_L = df_cell[k,:LB]#[row] 
                    c_U = df_cell[k,:UB]#[row] 
                    c = (c_U + c_L)/2
                    M = c_U - c_L
#                     println("cU = ", c_U)
#                     println("cL = ", c_L)
                    c_g = c + d.*x_now
                    y, gx, SP = gx_bound(c, c_g, edge)
                    df_cell[k,:g] = gx
                    df_cell[k,:Y] = y
#                     println("y = ", y)
                    println("Cell ",k,": ", findall(y.>0), " = ", gx)
#                     println("Shortest Path = ", findall(y.>0))
                    if  α_now - (z_now[k] + gx) > epsilon
                        
                        con_num = con_num + 1 
                        push!(df_constraints, (con_num, k, y, SP))
#                         println("epsilon-Condition not satisfied. Include y in path set for cell", k)
                        constr[con_num] = @constraint(m, α <= sum(d[i]*y[i]*x[i] for i = 1:Len) + SP + z[k])
#                         println(constr[con_num])
                        if α_now - (z_now[k] + gx) > delta1   
#                             println("O1 not satisfied: Add constraint for cell ", k)
                            O1Flag = false
                        end
                    end
                end
            end
            
            if O1Flag == true 
#                 O2Flag = true
#                 O3Flag = true
#                 partitionk = false
                
                #Verify against OC2
                for k in K_bar
#                     println("O2, ", k)
#                     row = findall(df_cell[!,:CELL] .== k)[1]
#                     c_L = df_cell[!,:LB][row] 
#                     c_U = df_cell[!,:UB][row]
                    c_L = df_cell[k,:LB]
                    c_U = df_cell[k,:UB]
                    M = c_U - c_L
                    c = (c_U + c_L)/2
                    c_g_L = c_L + d.*x_now
                    
                    #Solving for g(\hat{x}, c^{L,k})
                    yL, gL, SPL = gx_bound(c, c_g_L, edge)
                    df_cell[k,:Y_Lk] = yL
                    df_cell[k,:gL] = gL
#                     println("Cell ", k)
#                     println("α_now ", α_now)
#                     println("z_now[k]", z_now[k])
#                     println("gL = ", gL)
                    #Verify against OC2
                    if α_now - z_now[k] - gL <= delta2 #&& df_cell[k,:PROB] < β
#                         setdiff!(K_bar,k)
                        push!(K_removed,k)
                    else
#                         println("O2 not satisfied")
                        y_h, hx = hx_bound(c_L, c_U, d, x_now)
                        p_k = df_cell[k,:PROB]
                        df_cell[k,:h] = hx
                        gx = df_cell[k,:g] 
                        println("Cell ", k, ": gx = ", gx, "; hx = ", hx, " Y = ", findall(df_cell[k,:Y].>0))
                        
                        #Verify against OC3
                        if gx - hx <= delta2
                            push!(K_removed,k)
#                             setdiff!(K_bar,k) #Remove k from K_bar if satisfy OC3
                        else
                            
#                             println("Before calling Partition")
                            newCell = newCell + 1 #nrow(df_cell)+1
#                             newCell = maximum(df_cell[!,:CELL])+1
#                             K_partition, Δ, arc_split, yL, yU, gL, gU, SP_L, SP_U, O3Flag = Partition(df_cell, K_partition, newCell, k, p_k, gx, c_L, c_U, M, Y_k, O2Flag)
                            
                            Δ, arc_split, yL, yU, gL, gU, SP_L, SP_U = Partition(x_now, newCell, k, p_k, c_L, c_U, M, df_cell[k,:Y],yL)

                            #Updating constraints in MP:
#                             constr_of_k = findall(df_constraints[!,:CELL].== k) #Current constraints of k
                            
                            #Query rows from df_constraints that satisfy:
                            #(a) CELL column = k, and
                            #(b) Y column uses path that has arc # =  arc_split 
                            df_temp_k = df_constraints |> 
                            @filter(_.CELL == k)|> DataFrame
#                             @filter(_.Y[arc_split] == 1) |> DataFrame
                            
                            add_yL = true #Turns false if path yL is in cell k's existing constraints
                            add_yU = true #Turns false if path yU is in cell k's existing constraints
                            existingPath = false
                            
                            #LOOP THRU ALL ROWS IN DF ASSOCIATED WITH k
                            #Update RHS of affected constraints in k
                            #Copy each constraint in k to |K|+1
#                             println("Update constraints")
                            for dfRow in eachrow(df_temp_k) #constr_of_k 
#                                 Y_k = df_constraints[!,:Y][i] #Get path y in that constraint
                                Y_k = dfRow.Y
#                                 newCell_RHS = df_constraints[!,:SP][i] #RHS for newCell if arc_split is not in Y_k
                                if Y_k == yL #Found yL in the current P-set
                                    add_yL = false
                                end
                                if Y_k == yU #Found yU in the current P-set
                                    add_yU = false
                                end
                                
                                newCell_RHS = dfRow.SP
                                #ONLY UPDATE RHS IF ARC SPLIT IS ON THAT PATH  
                                if Y_k[arc_split] == 1
#                                 if Y_k[arc_split] == 1 #in findall(Y_k.==1)
#                                     k_new_RHS = df_constraints[!,:SP][i] - Δ #gap #FIND NEW RHS OF CELL k
                                    newCell_RHS = newCell_RHS + Δ #gap #FIND NEW RHS OF NEWCELL
                                    dfRow.SP = dfRow.SP - Δ
                                    conRef = dfRow.NUM
#                                     df_constraints[!,:SP][i] = k_new_RHS #FIX DF INFORMATION OF CELL k
                                    set_normalized_rhs(constr[conRef], dfRow.SP)
#                                     set_normalized_rhs(constr[i], k_new_RHS) #FIX CONSTRAINT RHS OF CELL k
                                end
                                
                                #COPY PATH Y_k FROM k to NEWCELL = |K|+1
                                con_num = con_num + 1 
                                constr[con_num] = @constraint(m, α <= 
                                            sum(d[i]*x[i]*Y_k[i] for i = 1:Len) + newCell_RHS + z[newCell])
                                push!(df_constraints, (con_num,newCell, Y_k, newCell_RHS)) #ADD DF INFORMATION OF CELL |K|+1
                            end

                            if add_yL == true #yL is a new path not in P^k
                                con_num = con_num + 1
                                constr[con_num] = @constraint(m, α <= 
                                            sum(d[i]*x[i]*yL[i] for i = 1:Len) + SP_L + z[k])
#                                 println("New path for cell k ", constr[con_num])
                                push!(df_constraints, (con_num,k, yL, SP_L))
                            end
                            if add_yU == true #yU is a new path not in P^{|K|+1}
                                con_num = con_num + 1
                                constr[con_num] = @constraint(m, α <= 
                                            sum(d[i]*x[i]*yU[i] for i = 1:Len) + SP_U + z[newCell])
#                                 println("New path for cell |K|+1 ", constr[con_num])
                                push!(df_constraints, (con_num, newCell, yU, SP_U))
                            end
                        end
                    end
                end #END OF for k = 1:myLength
                
                p = df_cell.PROB #[!,:PROB]
                @objective(m, Max, α - (1/(1-β))*sum(p[i]*z[i] for i = 1:length(p)))
                setdiff!(K_bar,K_removed)
                K_bar = vcat(K_bar, K_newly_added)
                K_newly_added = []
                K_removed = []
            end #If O1Flag == true
#         println("K_bar = ", K_bar)
        end #If Feasible
    end #While K_partition is non-empty
   
    if terminate_cond == false
        #Convolution
        println("Begin Convolution")
#         println(df_cell)
        for i in eachrow(df_cell)
           println("Cell ", i.CELL, " : ", findall(i.Y .> 0)) 
        end
        df_cellPoly = convolveEachCell()
#         println(df_cellPoly)
        #FindCVaR
        nu_L = minimum(df_cell.gL)
        yU, nu_U, SPU = gx_bound(cU_orig, cU_orig + d.*x_now, edge)
        println("Begin FindCVaR")
        CVaR = FindCVaR(α_now, nu_L, nu_U, df_cellPoly)
        K_bar = collect(1:nrow(df_cell))
        t_con = @constraint(m, sum(x_now[i]*x[i] for i = 1:Len) <= b-1)
        println("ADD X-CONSTRAINT TO MP: ", t_con, "\n")
        if LB < CVaR
#             LB = copy(CVaR)
            LB = CVaR
            x_sol, α_sol, z_sol = x_now, α_now, z_now
#             x_sol = copy(x_now) 
#             α_sol = copy(α_now)
#             z_sol = copy(z_now)
            
        end
        println("wth Incumbent solution:")
        println("x = ", findall(x_sol.==1))
        println("α_sol = ", α_sol)
#         println("z_sol = ", z_sol)
    end
    elapsed = time() - start
     
#     end
#     time_lapse = toq()    
#     total_time = total_time + time_lapse
#         if total_time > 3600
#             total_time = ">3600"
# #             println("total time" > 3600")
#             break
#         end
end
# println(df_cell)
# outfile = "C:/Users/din/Documents/GitHub/IntCVaR/df_cell.csv"
# f = open(outfile,"w") 
# println(f, df_cell)
# close(f)

# for i in eachrow(df_cell)
#     Y = findall(i.Y .>0)
#     D = d.*x_sol
#     println(i.LB[Y] + D[Y])
#     println(i.UB[Y] + D[Y])
# end
# CSV.write("df_cell.csv",Matrix(df_cell))
println("Terminating... ")
total_time = time() - start

# @timeit to "CALC OPT GAP" opt_gap = sum(p[i]*(z[i] - h[i]) for i = 1:length(p))
##println(f,"\n\n\nFinal z = ", getvalue(z[1:length(p)]))
# println(f,"Final g = ", g)
# println(f,"Final h = ",h)
##println(f,"List of all cells = ", Cell_List)

println("\n\nInstance ", myFile)
println("1-β = ", 1-β)
println("delta2 = ", delta2)
println("Interdiction = ", findall(x_sol.==1))
println("LB = ", LB, "; UB = ", MP_obj)
println("α_sol = ",α_sol)
println("Overall runtime = ", total_time)
println("No. iterations = ", iter) #length(df[:CELL]))

# println(elapsed)
