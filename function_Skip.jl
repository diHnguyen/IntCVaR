function Skip(K_not)
#     println("Inside function")
#     global α
#     println(df_Elim)
    K_α = []
    K_α_not = []
    K_not = []
    temp_df = sort(df_Elim, [:gU,:gL])
#     δ3 = 1.0
#     println("α = ", α)
    println(temp_df)
    temp_p = 0
    
    i=0
    terminate_cond = false
    while temp_p < α
        i = i+1
        temp_p = temp_p + temp_df[i,:PROB]
    end
    if temp_p < 1
        max_g = temp_df[i,:gU]
#         println("max_g = ", max_g)
#         println("i = ", i)
        cell_i = 0
#         K_α = temp_df[!,:CELL][1:i] #Group of cells making up alpha
        if temp_p > α
            K_α = temp_df[!,:CELL][1:(i-1)]
            cell_i = temp_df[!,:CELL][i]
        else
            K_α = temp_df[!,:CELL][1:i]
        end 
        
        myRows = collect((i+1):nrow(temp_df))
        min_g = minimum(temp_df[!,:gL][myRows])
        
        temp_K = []
        if min_g - max_g > δ3
            temp_K = filter(row -> row.gU > min_g, temp_df)[!,:CELL]  
        
    #         println("temp_K = ", temp_K)
    #         println("max_gU = ", max_g)
    #         println("min_gL = ", min_g)
    #         println("i = ", i)
    #         min_g = minimum(temp_df[myRows,:gL])
    #         println("min_gU = ", min_g)
            temp_df = temp_df[myRows,:]
    #         println(temp_df)
    #         println(filter!(row -> row.gL >= max_g, temp_df))
    #         println(filter!(row -> row.gL < max_g, temp_df))
            K_not = filter(row -> row.gL >= max_g, temp_df)[!,:CELL]   
    #         println("K_not = ", K_not)
    #         K_not = temp_df[!,:CELL] #Group of cells not in alpha that has gL > g_max 
            K_α_not = filter(row -> row.gL < max_g, temp_df)[!,:CELL] #Group of cells not in alpha
    #         println("K_α_not = ", K_α_not)
        
            if length(K_α_not) == 0#length(K_not) == nrow(df_Elim) - i
                terminate_cond = true
            else
                K_α_not = vcat(K_α_not,temp_K)
            end
        end
    
    else
        K_α = df_Elim[!,:CELL]
        K_not = []
    end
    
#     if length(K_not) == nrow(df_Elim) - i
#     println(temp_df)
#     push!(K_not,temp_df[!,:CELL])
    
#     println("K_α_not updated = ",K_α_not)
#     println("K_not length = ", length(K_not))
#     println("K_not weight total = ", sum(p[k] for k in K_not))
    return K_not, K_α, K_α_not, terminate_cond
    
end