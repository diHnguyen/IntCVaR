function chooseMax(M)
    #Choose an arc that has the max range of uncertainty
    arc_split = findall(M .== maximum(M))[1]
    return arc_split
end

function chooseScore(M, split_Candidates,k)
    #Choose an arc that has the lowest Score from Tree Sensitivity 
    #Analysis when Max_Range_Unc/Min_Range_Unc <= 10
    global c_L
    global c_U
    global c
    global c_g
    global x_now
    global edge
    
    Delta_Final_Pos = ones(Len).*-1
    Delta_Final_Neg = ones(Len).*-1

    M_max = maximum(M)
    temp = findall(M .>0)
    M_min = M_max
    for e in temp
        if M[e] <= M_min
            M_min = M[e]
        end
    end
    if M_max/M_min > 10
#                             println("Pick Max")
        split_Candidates[k] = findall(M .== maximum(M))[1]
    end
    if M_max/M_min <= 10
#                             println("SA")
        #FIND TREE, PRED
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

        arcs_neg = findall(Delta_Final_Neg .>= 0)
        arcs_pos = findall(Delta_Final_Pos .>= 0)
        found_split_arc = false
        arcs_delta = []
        myMax = maximum(Delta_Final_Neg)

        arcs = findall(M.>0)                        
        if myMax < maximum(Delta_Final_Pos)
            myMax = maximum(Delta_Final_Pos)
        end

        myScore = myMax/M_min
        for i in arcs
#                                 println("i = ", i)
            if Delta_Final_Pos[i] >= 0 && Delta_Final_Pos[i]/M[i] < myScore
                myScore = Delta_Final_Pos[i]/M[i]
                split_Candidates[k] = i
                found_split_arc = true
            end
            if Delta_Final_Neg[i] >= 0 && Delta_Final_Neg[i]/M[i] < myScore
                myScore = Delta_Final_Neg[i]/M[i]
                split_Candidates[k] = i
                found_split_arc = true
            end
        end
        if found_split_arc == false
            split_Candidates[k] = findall(M .== maximum(M))[1]
        end
    end
    return split_Candidates[k]
end

function chooseWorst(M, yL, yW)
    #After solving (3) at cL, solve (3) at cW where cW is worst for
    #yL. Compare yL and yW to obtain a set of candidate arcs. Choose
    #arc with largest uncertainty range.
    candidates = findall(yL.!=yW)
    arc_split = findall(M .== maximum(M[candidates]))[1]
    return arc_split
end
