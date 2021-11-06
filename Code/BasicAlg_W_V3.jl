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
global β # = rand(1:999)/1000
# gurobi_env.setParam("LogToConsole", 0)

to = TimerOutput()

dataSet = "N50_HalfUncRange"
Ins = string(ARGS[1])
myFile = "/Users/din/Desktop/IntCVaR/Instances/New Instances/"*dataSet*"/"*dataSet*"_"*Ins*".jl"
include(myFile)

include("functionGbound.jl")
include("functionHbound.jl")
include("functionPartition_BasicAlg.jl")
include("functionConvolution.jl")
include("functionFindCVaR_V3.jl")


println("Ins ", dataSet,"_", Ins, ": ", β, " Running...", myRun)
global epsilon = 1e-4
# global delta3 = 1.0
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

#Obtain the x* for worst case:
global x_sol = []
yL, gL, SP_L = gx_bound(cL_orig, cL_orig, edge)

m0 = Model(() -> Gurobi.Optimizer(gurobi_env)) # If we want to add # in Gurobi, then we have to turn of 
# set_optimizer_attribute(m, "OutputFlag", 0)    #Gurobi's own Cuts 
@variable(m0, x0[1:Len], Bin)
@variable(m0, z0)
@constraint(m0, sum(x0[i] for i=1:Len) == b) #attack budget = 7
@constraint(m0, z0 <= SP_L + sum(yL[i]*x0[i]*d[i] for i=1:Len)  )
@objective(m0, Max, z0)
global last_MP_0 = 1e6
global MP_0 = 1e5
optimize!(m0)
MP_0 = JuMP.objective_value.(m0)
global stop_det = false
while stop_det == false
   global MP_0, x_sol, stop_det
   stop_det = true
   last_MP_0 = JuMP.objective_value.(m0)
   x_sol = JuMP.value.(x0)
   z0_sol = JuMP.value.(z0)
   cLx = cL_orig + d.*x_sol
   y0, gx0, SP0 = gx_bound(cL_orig, cLx, edge)
   #println("gx0 = ", gx0, "; Obj = ", last_MP_0)
   if  z0_sol - gx0 > epsilon
      stop_det = false
      con = @constraint(m0, z0 <= sum(d[i]*y0[i]*x0[i] for i = 1:Len) + SP0)
      #println("Add con")
   end
   #println(m0)
   optimize!(m0) 
   MP_0 = JuMP.objective_value.(m0)
   #println("MP_0 = ", MP_0,"; lastMP = ", last_MP_0)
   x_sol = JuMP.value.(x0)
end

println("x_w = ", findall(x_sol .>0))


#MAIN PROGRAM:
df_constraints = DataFrame(NUM = Int[], CELL = Int[], Y = Array[], SP = Float64[])
df_cell = DataFrame(CELL = Int[], Y = Array[], Y_Lk = Array[], g = Float64[], h = Float64[], gL = Float64[], LB = Array[], UB = Array[], PROB = Float64[])


MP_obj = 0.0

push!(df_cell, (1, yy,yy, SP_init, 0, 0, cL_orig, cU_orig, 1))
push!(df_constraints, (1, 1,yy,SP_init))
##println(f,"MASTER PROBLEM==========================================================================================")

zNum = 100000
cRefNum = 2000000
m = Model(() -> Gurobi.Optimizer(gurobi_env)) # If we want to add # in Gurobi, then we have to turn of 
# set_optimizer_attribute(m, "OutputFlag", 0)    #Gurobi's own Cuts 
# println("1")
@variable(m, x[1:Len], Bin)
@variable(m, α)
@variable(m, 1e6 >= z[1:zNum] >= 0)
# @constraintref constr[1:200000]
# @ConstrRef constr[1:200000]
for a = 1:Len
    @constraint(m, x[a] == x_sol[a])
end

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
global LB_w = 0 
global MP_obj = 1e6
global K_newly_added = []
global K_removed = []
start = time()
global terminate_cond = false
global nu_U = 0
global nu_L = 1e6
global numConv = 0 
while terminate_cond == false 
    global α, β, iter, total_time, K_bar, K_newly_added, K_removed, LB, MP_obj, con_num, newCell , nu_U, nu_L, LB_w, numConv
    global x_sol, z_sol, α_sol, last_x, x_now,α_now,z_now, terminate_cond

    while isempty(K_bar) == false #length(K_bar) > K
        iter = iter + 1
        
        optimize!(m) 
#        println("\nIter : ", iter," ; LB = ", LB)
#         println(m)
#         println("", df_cell)
#         K = vcat(K, K_newly_added)
        

        if termination_status(m) == MOI.OPTIMAL
            MP_obj = JuMP.objective_value.(m)
            x_now = JuMP.value.(x)
            α_now = JuMP.value.(α)
            α_sol = α_now
            z_now = JuMP.value.(z)
            #println("\nIter : ", iter," ; Obj = ", MP_obj)
            #println(iter, " : ", length(K_bar),"/", newCell)
            #println("x = ", findall(x_now.>0))
        end
        
        if termination_status(m) != MOI.OPTIMAL || MP_obj <= LB
            terminate_cond = true
            K_bar = []
        else
            O1Flag = true
#            if last_x != x_now
                K_bar = collect(1:newCell)
		nu_U = 0
		nu_L = 1e6
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
#                     println("Cell ",k,": ", findall(y.>0), " = ", gx)
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
#            end
            
            if O1Flag == true 
                
                #Verify against OC2
#                 println("K_bar = ", K_bar)
                for k in K_bar
                    c_L = df_cell[k,:LB]
                    c_U = df_cell[k,:UB]
                    M = c_U - c_L
                    c = (c_U + c_L)/2
                    c_g_L = c_L + d.*x_now
		            yK = df_cell[k,:Y]
                    phi_yK_L = sum(c_g_L[i]*yK[i] for i = 1:Len)
                    phi_yK_U = phi_yK_L + sum(M[i]*yK[i] for i = 1:Len)
                    if phi_yK_L < nu_L
                        nu_L = phi_yK_L
                    end
                    if phi_yK_U > nu_U
                        nu_U = phi_yK_U
                    end
                    
                    #Solving for g(\hat{x}, c^{L,k})
                    yL, gL, SPL = gx_bound(c, c_g_L, edge)
                    df_cell[k,:Y_Lk] = yL
                    df_cell[k,:gL] = gL

                    #Verify against OC2
                    if α_now <= gL 
                        push!(K_removed,k)
                    else
                        y_h, hx = hx_bound(c_L, c_U, d, x_now)
                        p_k = df_cell[k,:PROB]
                        df_cell[k,:h] = hx
                        gx = df_cell[k,:g] 
                        println(k, ". ", gx, " -- ", hx)
                        #Verify against OC3
                        if gx - hx <= delta2
                            push!(K_removed,k)
                        else
                            
#                             println("Before calling Partition")
                            newCell = newCell + 1 #nrow(df_cell)+1
                            
                            Δ, arc_split, yL, yU, gL, gU, SP_L, SP_U = Partition(x_now, newCell, k, p_k, c_L, c_U, M, df_cell[k,:Y])

                            #Updating constraints in MP:
                            
                            #Query rows from df_constraints that satisfy:
                            #(a) CELL column = k, and
                            #(b) Y column uses path that has arc # =  arc_split 
                            df_temp_k = df_constraints |> 
                            @filter(_.CELL == k)|> DataFrame
                            add_yL = true #Turns false if path yL is in cell k's existing constraints
                            add_yU = true #Turns false if path yU is in cell k's existing constraints
                            existingPath = false
                            
                            #LOOP THRU ALL ROWS IN DF ASSOCIATED WITH k
                            #Update RHS of affected constraints in k
                            #Copy each constraint in k to |K|+1

                            for dfRow in eachrow(df_temp_k) #constr_of_k 
                                Y_k = dfRow.Y
                                if Y_k == yL #Found yL in the current P-set
                                    add_yL = false
                                end
                                if Y_k == yU #Found yU in the current P-set
                                    add_yU = false
                                end
                                
                                newCell_RHS = dfRow.SP
                                
                                #ONLY UPDATE RHS IF ARC SPLIT IS ON THAT PATH  
                                if Y_k[arc_split] == 1
                                    conRef = dfRow.NUM
                                    newCell_RHS = newCell_RHS + Δ #gap #FIND NEW RHS OF NEWCELL
                                    dfRow.SP = dfRow.SP - Δ
                                    
                                    df_constraints[conRef,:SP] = dfRow.SP
                                    set_normalized_rhs(constr[conRef], dfRow.SP)
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
                                push!(df_constraints, (con_num,k, yL, SP_L))
                            end
                            if add_yU == true #yU is a new path not in P^{|K|+1}
                                con_num = con_num + 1
                                constr[con_num] = @constraint(m, α <= 
                                            sum(d[i]*x[i]*yU[i] for i = 1:Len) + SP_U + z[newCell])
                                push!(df_constraints, (con_num, newCell, yU, SP_U))
                            end
                        end
                    end
                end #END OF for k = 1:myLength
                
                p = df_cell.PROB #[!,:PROB]
                #@objective(m, Max, α - sum(p[i]*z[i] for i = 1:length(p)))
                @objective(m, Max, α - (1/(1-β))*sum(p[i]*z[i] for i = 1:length(p)) )
                setdiff!(K_bar,K_removed)
                K_bar = vcat(K_bar, K_newly_added)
                K_newly_added = []
                K_removed = []
            end #If O1Flag == true
        end #If Feasible
    end #While K_partition is non-empty
    println("Convolve x_now = ", findall(x_now.>0))
    if terminate_cond == false
	numConv = numConv+1
        terminate_cond = true
        #Convolution
        df_cellPoly = convolveEachCell()
        println("nu_L = ", nu_L)
        println("nu_U = ", nu_U)
        α_now = 0.5*(nu_L + nu_U)
        CVaR, weight = FindCVaR(α_now, nu_L, nu_U, df_cellPoly)
        K_bar = collect(1:nrow(df_cell))

        #@constraint(m, sum(x_now[i]*x[i] for i = 1:Len) <= b-1)
        #println("ADD X-CONSTRAINT TO MP: ", t_con, "\n")

        println("CVaR = ", CVaR)
        LB = CVaR
	LB_w = weight
        x_sol, z_sol = x_now, z_now
    end
end

total_time = time() - start

println("BasicAlgW_V3_"*dataSet, "; Ins ", Ins, "; β ",β,"; Time ", total_time, "; LB ", LB, "; α_sol ",α_sol,"; x_sol ", findall(x_sol.==1),"; Cells ", nrow(df_cell), "; Iter ", iter, "; W ", LB_w, "; Cuts ", numConv)

timesFile = open("/Users/din/Desktop/IntCVaR/OutputFile/BasicAlgW_V3_"*dataSet*".txt", "a")
println(timesFile, dataSet, "; Ins ", Ins, "; β ",β,"; Time ", total_time, "; LB ", LB, "; α_sol ",α_sol,"; x_sol ", findall(x_sol.==1),"; Cells ", nrow(df_cell), "; Iter ", iter, "; W ", LB_w, "; Cuts ", numConv)
close(timesFile)
println(LB_w + β)
println("\007")