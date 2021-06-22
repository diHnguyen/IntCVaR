#Ready for upload
using JuMP 
using Gurobi
using LightGraphs
using DataFrames
using CSV
using TimerOutputs
using Dates

myRun = Dates.format(now(), "HH:MM:SS")
global gurobi_env = Gurobi.Env()
# gurobi_env.setParam("LogToConsole", 0)

to = TimerOutput()
myFile = "./myData1.jl"
include(myFile)
include("PartitionRules.jl")
include("functionGbound.jl")
include("functionHbound.jl")
include("functionPartition_BasicAlg.jl")
include("functionConvolution.jl")

println("Running...")
global epsilon = 1e-3
global delta3 = 5
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
df_cell = DataFrame(CELL = Int[], Y = Array[], Y_Lk = Array[], g = Float64[], h = Float64[], gL = Float64[], LB = Array[], UB = Array[], PROB = Float64[])


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
global x_sol = []
global α_sol = 0
global z_sol = []
global x_now = []
global α_now = 0
global z_now = []
global last_x = zeros(Len)
global con_num = 1
global total_time = 0.0
global iter = 0
global K = Int64[1]
global K_bar = Int64[]
global LB = 0
global MP_obj = 1e6
start = time()
terminate_cond = false

while terminate_cond == false 
    global α, β, iter, total_time, K, K_removed, LB, MP_obj, con_num, min_gLk, max_gUk
    global x_sol, z_sol, α_sol, last_x, x_now,α_now,z_now
    println("lengthK = ", length(K))
    while length(K_bar) > K
        stopping_cond = 0
        iter = iter + 1
        optimize!(m) 
        println("\nIter : ", iter)
        
        if termination_status(m) == MOI.OPTIMAL
            MP_obj = JuMP.objective_value.(m)
            x_now = JuMP.value.(x)
            println("Obj = ", MP_obj)
            println("Interdiction ", findall(x_now.==1))
            α_now = JuMP.value.(α)
            last_x = x_now
            z_now = JuMP.value.(z)
        end
        
        if termination_status(m) != MOI.OPTIMAL || MP_obj <= LB
            terminate_cond = true
            println("MP not feasible... Terminating")
        else
            O1Flag = true
            if last_x != x_now
                K_bar = collect(1:myLength)
                last_x = x_now
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
                    if  α_now - (z_now[k] + gx) > epsilon
                        push!(df_constraints, (k, y, SP))
                        con_num = con_num + 1 
#                         println("epsilon-Condition not satisfied. Include y in path set for cell", k)
                        constr[con_num] = @constraint(m, α <= sum(d[i]*y[i]*x[i] for i = 1:Len) + SP + z[k])
                        println(constr[con_num])
                        if α_now - (z_now[k] + gx) > delta1   
                            println("O1 not satisfied: Add constraint for cell ", k)
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
                    row = findall(df_cell[!,:CELL] .== k)[1]
                    c_L = df_cell[!,:LB][row] 
                    c_U = df_cell[!,:UB][row]
                    c = (c_U + c_L)/2
                    M = c_U - c_L
                    c_g_L = c + d.*x_now
                    
                    #Solving for g(\hat{x}, c^{L,k})
                    yL, gL, SPL = gx_bound(c, c_g_L, edge)
                end
                
                #STOPPED HERE LAST ON 06/22
                
                
                
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
#                             df_cell[k,:gL] = gL
                            if α_now - gL > δ3
                                println("+ Must partition")
                                O3Flag = false
                            end
                        else

                            gL = sum(c_L[i]*Y_k[i] for i = 1:Len)
                            gU = sum(c_U[i]*Y_k[i] for i = 1:Len)
                            if gL < min_gLk
                                min_gLk = gL + 0.0
                            end
                            if gU > max_gUk
                               max_gUk = gU + 0.0
                            end
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
                if !(length(p) in K)
                    push!(K, length(p))
                end
            end #If O1Flag == true
            
        println("K = ", K)
        println("K_removed = ", K_removed)
#         filter!(x-> !(x in K_removed), K)
#         K_removed = []
#         println("K = ", K)
#         println("K_removed = ", K_removed)
        end #If Feasible
#         println("Set of newly added cells: ", K_newCells)
#         K_partition = vcat(K_partition, K_newCells)
    end #While K_partition is non-empty
    println(min_gLk,", ",max_gUk)
    println("df_cell = ")
    println(df_cell)
    if length(K) == length(K_removed)#isempty(K) == true
        @constraint(m, sum(x_now[i]*x[i] for i=1:Len) <= b-1)
        println("Check beta-CVar")
        df_cellPoly = DataFrame(CELL = Int64[], df = Any[], NUMPOLY = Int64[], DETSHIFT = Float64[])
        K_removed = [] #copy(K)
#         break
        println("K = ", K)
        for k in K# 1:3 #length(df_cell)
            cL = df_cell[k,:LB]
            cU = df_cell[k,:UB]
            y = df_cell[k,:Y]
            Len = length(cL)
            println("\nCell ", k)
            unc_arcs = findall((y.>0) .& (cL.<cU))
            det_arcs = findall((y.>0) .& (cL.==cU))
            println("unc_arcs = ", unc_arcs)
            println("det_arcs = ", det_arcs)
            
            if isempty(det_arcs) == false
                det_Shift = sum(cL[i] for i in det_arcs)
            else
                det_Shift = 0
            end
        #     det_arcs = findall((Y.>0) .& (cL.==cU))
        #     filter!(x-> !(x in unc_arcs), det_arcs)
        #     det_Shift = 0
            if length(unc_arcs) > 0
                
                df, numPoly = Convolve(cL, cU, Len, unc_arcs)

            #     df[:,:leftShift] = df[:,:leftShift] .+ det_Shift
                df.w = zeros(numPoly)
                println(df)
                push!(df_cellPoly, (k,df,numPoly, det_Shift))
            else
                push!(df_cellPoly, (k,DataFrame(),-1, det_Shift))
            end
        end
        println(df_cell[:,:PROB])
        println(df_cellPoly)
#         break
        temp_CVaR = FindCVaR(β,α_now,min_gLk,max_gUk,df_cellPoly, df_cell[:,:PROB])
        if LB < temp_CVaR 
            LB = temp_CVaR    
            x_sol = x_now
            α_sol = α_now
            z_sol = z_now
        end
    end
    if iter > 6
        println("Break early...")
        terminate_cond = true
#             K = []

        K_removed = copy(K)
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
