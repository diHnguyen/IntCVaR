include("./functions_sharedTreePath.jl")
function what_arc(T, y, e)
    arc_type = -1 #type 0 = basic-SP, 1 = basic-nonSP, 2 = non-basic
    if T[e] == 1
        if y[e] == 1
            arc_type = 0
        else
            arc_type = 1
        end
    else
        arc_type = 2
    end
end

function find_affected_node(T, pred, e, edge) #NEED TO CONSIDER CYCLES
    S = Int[]
    N = Int[]
    push!(S, edge[e,2])
    push!(N, edge[e,2])
    visit = collect(1:last_node)
    l = length(S)
    visit[edge[e,2]] = 0
    while l > 0 
        node = S[1]
        for i = 1:last_node 
            if pred[i] == node && visit[i] != 0
                push!(S, i)
                push!(N, i)
                visit[i] = 0
            end
            
        end
        S[1], S[length(S)] = S[length(S)], S[1]
        pop!(S)
        l = length(S)
    end
    return N
end

function SA_basic(N, nonbasic, nonbasic_cost, label_temp, temp)
    delta = []
    delta_ind = []
    δ_ind = 0
    val_temp = 0
    if length(N) > 0
        ##println("nonbasic = ", nonbasic)
        for j in N
            tail = findall(nonbasic[:,1] .== j)
            if length(tail) > 0
                ##println("tail = ", tail)
                for i in tail
                    ##println("Arc #", i)
                    ##println("Arc name = ", nonbasic[i,:])
                    condd = nonbasic[i,2] in N
                    if condd == false
                        ##println("FALSE THEN")
                        val_temp = label_temp[nonbasic[i,1]] + nonbasic_cost[i] - label_temp[nonbasic[i,2]]
                        push!(delta, val_temp)
                        push!(delta_ind, i)
                    end
                end
            end
        end
    end
    
    ##println("delta = ", delta)
    if length(delta) > 0
        ##println("DELTA > 0")
        for i = 1:length(delta)
            if i == 1
                δ_min = delta[1]
                δ_ind = delta_ind[1]
            end
            if delta[i] < δ_min
                δ_min = delta[i]
                δ_ind = delta_ind[i]
            end
        end
        temp = temp - δ_min
    else
        δ_min = -1.0
    end

    return δ_min, δ_ind, temp
end

function find_SA_Delta(myNode, Delta_Final_Pos, Delta_Final_Neg, sign, basic, i_basic, nonbasic, i_nonbasic)
    myNonTreeArcList = []
    myTreeArcList = []
    Delta_List = []
    delta = maximum(Delta_Final_Neg)
    entering_arc = 0
    found = false
#     println("myNode = ", myNode)
#     println(nonbasic
    if sign == "POS"
        myNonTreeArcList = findall(nonbasic[:,2].== myNode)
        Delta_List = Delta_Final_Pos
    end
    if sign == "NEG"
        myNonTreeArcList = findall(nonbasic[:,1].== myNode)
        Delta_List = Delta_Final_Neg
    end
    
    #NONBASIC ARCS
    if length(myNonTreeArcList) > 0
#         println("Nonbasic")
#         println("myNonTreeArcList = ", myNonTreeArcList)
        
#         e = i_nonbasic[myNonTreeArcList[1]]
#         println("Arc ", e)
#         delta = Delta_Final_Neg[e]
#         println("New delta = ", delta)
        for i = 1 : length(myNonTreeArcList)
            e = i_nonbasic[myNonTreeArcList[i]]
#             println("Arc ", e, "Delta ", Delta_Final_Neg[e])
            if delta >= Delta_Final_Neg[e] && Delta_Final_Neg[e] >= 0
                delta = Delta_Final_Neg[e]
                entering_arc = e
#                 println("New delta = ", delta)
            end
        end
    end
    
    
    #BASIC ARCS
    myTreeArcList = findall(basic[:,1].== myNode)
#     println("basic = ", basic)
#     println("myTreeArcList = ", myTreeArcList)
#     println("l = ",length(myTreeArcList))
    
#     println("myTreeArcList = ", myTreeArcList)
    if length(myTreeArcList) > 0
#         println("Basic")
        if delta == -1
#             e = i_basic[myTreeArcList[1]]
            delta = maximum(Delta_List)
#             println("New delta = ", delta)
        end
        for i = 1 : length(myTreeArcList)
            e = i_basic[myTreeArcList[i]]
#             println("Arc ", e, "Delta ", Delta_List[e])
            if delta >= Delta_List[e] && Delta_List[e] >= 0
                delta = Delta_List[e]
                entering_arc = e
#                 println("New delta = ", delta)
            end
        end
    end
#     println("entering arc inside function: ", entering_arc)
    return delta, entering_arc
end

function findNonBasicArc_Loop(label, T, c_g, e, path, y, pred, basic, i_basic)
    label_temp = copy(label)
    Δ = 0
    δ = []
    δ_ind = 0
    T_temp = copy(T)
    pred_temp = copy(pred)
    c_temp = copy(c_g)
    edge_num = copy(e)
    N = Int[]
    cur_edge = Int[]
    path_temp = copy(path)
    
    new_SP = false
    
    
    u = edge[e,1]
    v = edge[e,2]
#     t = pred[v]
#     e_temp = [t v]
    ###println(f,"(t v) = (",t, " ", v, ")")
    i = 0
    for temp_i in i_basic
        if edge[temp_i, 1] == pred[v] && edge[temp_i,2] == v
            i = temp_i
            break
        end
    end
#     B = find(all(basic .== e_temp,2)) #find a tree edge to be replaced by nonbasic e
#     println("basic = ", basic)
#     println("e_temp = ", e_temp)
#     i = B[1]
#     println("Entering arc = ", e )
#     println("Exiting arc = ", i)
#     i = i_basic[i]
    val = label_temp[u] - label_temp[v] + c_temp[e]
    push!(δ, val)
    ###println(f,"δ = ", δ)
    δ_min = minimum(δ)
    c_temp[e] = c_temp[e] - δ_min
    Δ = Δ + δ_min
    pred_temp[v] = u
    ###println(f, "Tree before replace: ", T_temp)
    T_temp[e] = 1
    T_temp[i] = 0
    delta = -0.5
    while new_SP == false
        basic, basic_cost, i_basic, nonbasic, nonbasic_cost, i_nonbasic = classify_edges_costs(T_temp,c_temp)
        #loop = loop+1
        ###println("Loop ", loop)
        ##println( "Current T = ", T_temp)
        ##println( "Current y = ", y)
        ##println("Current label = ", label_temp)
        arc_type = what_arc(T_temp, y, e)
#         println("Double check: Arc ", e, " is now type ", arc_type)
        if arc_type == 0 #THIS STATEMENT IS WORKING
            new_SP = true
            delta = Δ
        end
        if arc_type == 1
            N = find_affected_node(T_temp, pred_temp, e, edge)
            ##println("N = ", N)
            δ_min, δ_ind, temp = SA_basic(N, nonbasic, nonbasic_cost, label_temp, c_temp[e])
            ##println("Δ = ", Δ)
            ##println("δ_min = ", δ_min)
            if δ_min < 0
                delta = -0.5 #δ_min
                new_SP = true
            else
                ##println("nonbasic arcs: ", nonbasic)
                ##println("index of arc moving into the basis: ", δ_ind)
                u = nonbasic[δ_ind, 1]
                v = nonbasic[δ_ind, 2]
                t = pred_temp[v]
                e_temp = [t v]
                ##println("Basic arc becoming nonbasic: ", e_temp)
                ##println("Nonbasic arc becoming basic: ", [u v])
                B = find(all(basic .== e_temp,2)) #find a tree edge to be replaced by nonbasic e
                i = B[1]
                i = i_basic[i]
                pred_temp[v] = nonbasic[δ_ind, 1]
                c_temp[e] = temp
                T_temp[i_nonbasic[δ_ind]] = 1
#                 println("Entering arc ", i_nonbasic[δ_ind])
#                 println("Exiting arc ", i)
                T_temp[i] = 0 ################
                Δ = Δ + δ_min
                delta = Δ
#                 println("current delta = ", delta)
                ##println("Pred = ", pred_temp)
            end
            
        end 
        #CHECK FOR LOOP:


        if y.*T_temp != y
            new_SP = true
        else
            loop = false
            if pred_temp[origin] != 0
                loop = true
            else
                v = edge[e,2]
                visit = collect(1:last_node)
                ##println("visit = ", visit)
                while v != origin && visit[v] != 0
                    ##println("v = ", v)
                    ##println("pred[v] = ", pred_temp[v] )
                    visit[v] = 0
                    v = pred_temp[v]

                    ##println("visit = ", visit)

                end

                if v != origin
                    loop = true
                    
                    ##println("HERE")
                end
            end
            if loop == true
                new_SP = true

            end
        end
        ##println("new_SP = ", new_SP)

    end

    return delta
end