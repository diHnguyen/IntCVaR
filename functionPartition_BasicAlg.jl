# using JuMP
function Partition(df_cell, K_removed, newCell, k, p, gx, cL, cU, M, y, O2Flag)
#     println(O2Flag)
    global Len
    if O2Flag == true
        My = M.*y
        arc_split = findall(My.==maximum(My))[1]
        
        yL = y
        yU = y
        
        #We need to think about h-value for the new cells. Do we need to update h?
        #Most likely yes -- Remove the update on h-values below when figure out what to do with it. 
        gap =  df_cell[k,:g] - df_cell[k,:h]
        #====#
        Δ = M[arc_split]/4
        SP = sum(y[i]*(cL[i]+cU[i])/2 for i = 1:Len)
        SP_L = SP - Δ
        SP_U = SP + Δ
        gL = gx - Δ
        gU = gx + Δ
        
        M_ = zeros(Len)
        M_[arc_split] = M[arc_split]/2
        
#         cMid = copy(cL)
#         cMid[arc_split] = (cL[arc_split] + cU[arc_split])/2
        push!(df_cell, (newCell, yU, gU, gU, cL+M_, cU, p/2))
        
        df_cell[k,:UB] = cU - M_
        df_cell[k,:g] = gL
        df_cell[k,:h] = gL
        df_cell[k,:PROB] = p/2
        
        if gU - α_now > δ3 
            push!(K_removed, newCell)
        end
        if α_now - gL > δ3 
            push!(K_removed, k)
        end
    else
        arc_split = findall(M.==maximum(M))[1]
        M_ = zeros(Len)
        M_[arc_split] = M[arc_split]/2
        Δ = M[arc_split]/4
#         println("cU = ", cU)
#         println("c+ = ", cL+M_)
#         println("c- = ", cU-M_)
#         println("cL = ", cL)
#         println("arc_split = ", arc_split)
        
        cL_avg = (cL+cU-M_)/2
        cU_avg = (cL+cU+M_)/2
        yL, gL, SP_L = gx_bound(cL_avg, cL_avg+d.*x_now, edge)
        yU, gU, SP_U = gx_bound(cU_avg, cU_avg+d.*x_now, edge)
        
#         println("g revised k = ", gL)
#         println("g  = ", gU)
        
        push!(df_cell, (newCell, yU, gU, 0, cL+M_, cU, p/2))
        df_cell[k,:UB] = cU-M_
        df_cell[k,:g] = gL
        df_cell[k,:h] = 0
        df_cell[k,:PROB] = p/2
#         println("Done.")
    end
#     push!(K_partition, newCell)
    return K_removed,Δ, arc_split, yL, yU, gL, gU, SP_L, SP_U, false
end
