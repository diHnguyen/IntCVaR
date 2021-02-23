function classify_edges_costs(T,c)
    basic = Array{Int64}(undef, 0, 2)
    basic_cost = []
    i_basic = []
    nonbasic = Array{Int64}(undef, 0, 2)
    nonbasic_cost = []
    i_nonbasic = []
    
    for e = 1:Len
        #arc_type = what_arc(T,y,e)
        if T[e] != 1
            nonbasic = vcat(nonbasic, edge[e:e,:])
#             nonbasic = [nonbasic; edge[e:e,:]]
#             println(nonbasic)
            push!(i_nonbasic, e)
            push!(nonbasic_cost, c[e])
        else 
            basic = vcat(basic, edge[e:e,:])
#             basic = [basic; edge[e:e,:]]
            push!(i_basic, e)
            push!(basic_cost, c[e])
        end
    end
    #We need to index back to the original edge values for reference purposes
    return basic, basic_cost, i_basic, nonbasic, nonbasic_cost, i_nonbasic
end


function getSortedBasicArcList(pred, basic, i_basic)
#     succ = pred.*0
    sorted_basic = copy(basic)
    sorted_i_basic = zeros(Int64,length(i_basic))
#     println("pred = ", pred)
#     println("basic = ", basic)
#     println("i_basic = ", i_basic)
    
    index_sorting = length(i_basic)
#     println(index_sorting)
    S = [origin]
#     println("S = ",S)
    while isempty(S) == false
#         println("S = ",S)
        tailNode = S[1]
        deleteat!(S, 1)  
        for i = 1:length(pred)
            #origin vs destination
            if i != tailNode && pred[i] == tailNode
                u, v = tailNode, i
#                 println("(",u , ", ", v,")")
#                 myArc = transpose([u, v])
                temp_index = 0 #find(all(basic .== myArc,2))
#                 println("basic = ", basic)
                for j in i_basic
#                     println(edge[j,:], ": ", u, ", ", v)
                    if edge[j,1] == u && edge[j,2] == v #e_temp
                        temp_index = j
                        break
                    end
                end
#                 println(temp_index)
                if temp_index != 0#length(temp_index) > 0
#                     sorted_i_basic[index_sorting] = i_basic[temp_index[1]]
                    sorted_i_basic[index_sorting] = temp_index
                    index_sorting = index_sorting - 1
                    push!(S,i)
                end
#                 println("sorted_i_basic = ", sorted_i_basic)
            end
        end   
#         println("End of for loop")
    end
    return sorted_i_basic
#     println("sorted_i_basic = ", sorted_i_basic)
end
