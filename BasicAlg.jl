#Ready for upload
using JuMP 
using Gurobi
using LightGraphs
using DataFrames
using CSV
using TimerOutputs
using Dates

myRun = Dates.format(now(), "HH:MM:SS")
const gurobi_env = Gurobi.Env()
# gurobi_env.setParam("LogToConsole", 0)

to = TimerOutput()
# myFile = "./SmallExample.jl"
# myFile = "./testing.jl"
myFile = "./myData1.jl"
# myInstance = 7
# myFile = "./Instances/testInstance_"*string(myInstance)*".jl"
include(myFile)
include("PartitionRules.jl")
include("functionGbound.jl")
include("functionHbound.jl")
include("functionPartition_BasicAlg.jl")
include("functionConvolution.jl")

println("Running...")
epsilon = 1e-3
δ3 = 5
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
df_constraints = DataFrame(CELL = Int[], Y = Array[], SP = Float64[])
df_cell = DataFrame(CELL = Int[], Y = Array[], g = Float64[], h = Float64[], LB = Array[], UB = Array[], PROB = Float64[])
global last_x = zeros(Len)

MP_obj = 0.0

push!(df_cell, (1, yy, SP_init, 0, cL_orig, cU_orig, 1))
push!(df_constraints, (1,yy,SP_init))
##println(f,"MASTER PROBLEM==========================================================================================")

zNum = 100
cRefNum = 200
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
@constraint(m, sum(x[i] for i=1:Len) <= b) #attack budget = 2
constr[1] = @constraint(m, α <= SP_init + sum(yy[i]*x[i]*d[i] for i=1:Len) + z[1]  )
@objective(m, Max, α - (1/(1-β))*sum(p[i]*z[i] for i = 1:length(p)) )#w - sum(s[k] for k=1:length(s))/length(s) )


global con_num = 1
global total_time = 0.0
global iter = 0
global K = Int64[1]
global K_removed = Int64[]
start = time()
terminate_cond = false

while terminate_cond == false 
    println("lengthK = ", length(K))
    while length(K) > 0 #original : length(K_partition)
        stopping_cond = 0
        
        global iter = iter + 1
        
        global con_num
        global total_time
        global K
        global K_removed
        global last_x
        global x_now
        if iter > 10
            println("Break early...")
            terminate_cond = true
            K = []
            break
        end
#         println(m)
#         set_optimizer_attribute(m, "OutputFlag", 0)
#         println(m)
#         for i = 1:con_num
#            println(constr[i]) 
#         end
#         println(df_constraints)
        
#         println(df_cell)
        
        optimize!(m) 
        K_newCells = []
        println("\nIter : ", iter)
        
        if termination_status(m) != MOI.OPTIMAL
            terminate_cond = true
            println("MP not feasible... Terminating")
        else
            MP = JuMP.objective_value.(m)
            MP_obj = copy(MP)
            x_now = JuMP.value.(x) #+ 0.0 

            
            println("Obj = ", MP_obj)
            println("Interdiction ", findall(x_now.==1))
            #println("# Cells = ", length(p))
            α_now = JuMP.value.(α)
            last_x = x_now
            z_now = JuMP.value.(z)
            
        
            O1Flag = true
            if last_x != x_now
                K_partition = collect(1:myLength)
                last_x = x_now
            end
#             println("K_partition - ", K_partition)
            for k in K #_partition  
                row = findall(df_cell[!,:CELL] .== k)[1]
                c_L = df_cell[!,:LB][row] 
                c_U = df_cell[!,:UB][row] 
                c = (c_U + c_L)/2
                M = c_U - c_L
                println("cU = ", c_U)
                println("cL = ", c_L)
                c_g = c + d.*x_now
                y, gx, SP = gx_bound(c, c_g, edge)
                df_cell[k,:g] = gx
                df_cell[k,:Y] = y
                println("y = ", y)
                if  α_now - (z_now[k] + gx) >= epsilon
                    push!(df_constraints, (k, y, SP))
                    con_num = con_num + 1 
                    println("O1 not satisfied: Add constraint for cell ", k)
                    constr[con_num] = @constraint(m, α <= sum(d[i]*y[i]*x[i] for i = 1:Len) + SP + z[k])
                    println(constr[con_num])
                    if α_now - (z_now[k] + gx) > delta1     
                        O1Flag = false
                    end
                end
            end
            
            
            if O1Flag == true 
                O2Flag = true
                O3Flag = true
                partitionk = false
#                 K_partition_Len = length(K_partition)
                
                for k in K 
#_partition# 1:myLength #length(p)# k = 1: cell_num #length(p) 
#                     k = K_partition[ii]
#                     println("\n",k,". K_partition ", K_partition)
                    
                    println("K_removed = ", K_removed)
                    if (k in K_removed) == false
                        println("Partition? Cell ", k)

                        row = findall(df_cell[!,:CELL] .== k)[1]
                        c_L = df_cell[!,:LB][row] 
                        c_U = df_cell[!,:UB][row]
                        c = (c_U + c_L)/2
                        M = c_U - c_L
                        c_g = c + d.*x_now
                        Y_k = df_cell[k,:Y]
                        println("Y_k = ", Y_k)
                        p_k = df_cell[k,:PROB]
                        ##println(f,"==================SUBPROBLEM 2: H(X), CURRENT CELL = ", k)
                        y_h, hx = hx_bound(c_L, c_U, d, x_now)

                        df_cell[k,:h] = hx
                        gx = df_cell[k,:g] 
                        println("gx = ", gx, "; hx = ", hx)
                        if gx - hx > delta2
                            println("O2 not satisfied")
                            O2Flag = false
                            y, gL = gx_bound(c_L, c_L+d.*x_now, edge)
                            if α_now - gL > δ3
                                println("+ Must partition")
                                O3Flag = false
                            end
                        else

                            gL = sum(c_L[i]*Y_k[i] for i = 1:Len)
                            gU = sum(c_U[i]*Y_k[i] for i = 1:Len)
                            println("O2 is satisfied: Checking gL , gU")
                            println("gU = ", gU, "; ", "α_now = ", α_now,"; gL = ", gL)
                            if gU - α_now > δ3 && α_now - gL > δ3 
                                println("O3 not satisfied.")
                                O3Flag = false
                            end
                        end

                        if O3Flag == false
                            println("Before calling Partition")
                            newCell = maximum(df_cell[!,:CELL])+1

                            push!(K_newCells, newCell)

#                             K_partition, Δ, arc_split, yL, yU, gL, gU, SP_L, SP_U, O3Flag = Partition(df_cell, K_partition, newCell, k, p_k, gx, c_L, c_U, M, Y_k, O2Flag)
                            K_removed, Δ, arc_split, yL, yU, gL, gU, SP_L, SP_U, O3Flag = Partition(df_cell, K_removed, newCell, k, p_k, gx, c_L, c_U, M, Y_k, O2Flag)

        #                     println("Done partition")
                            constr_of_k = findall(df_constraints[!,:CELL].== k)

                            add_yL = true
                            add_yU = true

                            if O2Flag == true
                                add_yL = false
                                add_yU = false
                            end
                            #LOOP THRU ALL ROWS IN DF ASSOCIATED WITH k
                            #Update RHS of affected constraints in k
                            #Copy each constraint in k to |K|+1
                            println("Update constraints")
                            for i in constr_of_k 
                                Y_k = df_constraints[!,:Y][i] #Get path y in that constraint
                                newCell_RHS = df_constraints[!,:SP][i] #RHS for newCell if arc_split is not in Y_k
                                if Y_k == yL
                                    add_yL = false
                                end
                                if Y_k == yU
                                    add_yU = false
                                end
                                #ONLY UPDATE RHS IF ARC SPLIT IS ON THAT PATH           
                                if arc_split in findall(Y_k.==1)

                                    k_new_RHS = df_constraints[!,:SP][i] - Δ #gap #UPDATE RHS OF CELL k
                                    newCell_RHS = df_constraints[!,:SP][i] + Δ #gap #UPDATE RHS OF NEWCELL

                                    df_constraints[!,:SP][i] = k_new_RHS #FIX DF INFORMATION
                                    set_normalized_rhs(constr[i], k_new_RHS)
                                end
                                #COPY THE CONSTRAINT FROM k to NEWCELL
                                con_num = con_num + 1
                                constr[con_num] = @constraint(m, α <= 
                                            sum(d[i]*x[i]*Y_k[i] for i = 1:Len) + newCell_RHS + z[newCell])
            #                                 push!(df, (newCell, ARCS_USED_i, c_Mid_newCell, c_U, p_temp, newCell_RHS))
                                push!(df_constraints, (newCell, Y_k, newCell_RHS) )
                            end

                            if add_yL == true #yL is a new path not in P^k
                                con_num = con_num + 1
                                constr[con_num] = @constraint(m, α <= 
                                            sum(d[i]*x[i]*yL[i] for i = 1:Len) + SP_L + z[newCell])
                                println("New path for cell k ", constr[con_num])
                                push!(df_constraints, (k, yL, SP_L))
                            end
                            if add_yU == true #yU is a new path not in P^{|K|+1}
                                con_num = con_num + 1
                                constr[con_num] = @constraint(m, α <= 
                                            sum(d[i]*x[i]*yU[i] for i = 1:Len) + SP_U + z[newCell])
                                println("New path for cell |K|+1 ", constr[con_num])
                                push!(df_constraints, (newCell, yU, SP_U))
                            end
                        end
                        if O3Flag == true
                            println("O3 is satisfied")
                            push!(K_removed, k)
#                             println("K_partition = ", K_partition)
    #                         deleteat!(K, K.==2)
                            #Using filter! is much faster than deleteat! 
#                             filter!(x->x != k,K_partition)
    #                         println("k in K_par = ", (k in K_partition))
#                             deleteat!(K_partition, K_partition .== k);
#                             println("K_partition updated ", K_partition)
                        end
                    end
                end #END OF for k = 1:myLength
                p = df_cell[!,:PROB]
                @objective(m, Max, α - (1/(1-β))*sum(p[i]*z[i] for i = 1:length(p)))
                println("Update Obj Function")
            end #If O1Flag == true
            
        println("K = ", K)
        println("K_removed = ", K_removed)
        filter!(x-> !(x in K_removed), K)
        K_removed = []
        println("K = ", K)
        println("K_removed = ", K_removed)
        end #If Feasible
#         println("Set of newly added cells: ", K_newCells)
#         K_partition = vcat(K_partition, K_newCells)
    end #While K_partition is non-empty
    
    if isempty(K) == true
        println("Check beta-CVar")
        
        terminating_cond = true
        break
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

println("Terminating... ")
total_time = time() - start

@timeit to "CALC OPT GAP" opt_gap = sum(p[i]*(z[i] - h[i]) for i = 1:length(p))
##println(f,"\n\n\nFinal z = ", getvalue(z[1:length(p)]))
# println(f,"Final g = ", g)
# println(f,"Final h = ",h)
##println(f,"List of all cells = ", Cell_List)

println("\n\nInstance ", myFile)
println("1-β = ", 1-β)
println("delta2 = ", delta2)
println("Interdiction = ", findall(last_x.==1))
println("MP_obj = ", MP_obj)
println("K_not length = ", length(K_not), " out of ", length(p))
if isempty(K_not)
    println("K_not weight total = ", 0)
else
    println("K_not weight total = ", sum(p[k] for k in K_not))
end
println("Overall runtime = ", total_time)
println("No. iterations = ", iter) #length(df[:CELL]))

println(elapsed)
