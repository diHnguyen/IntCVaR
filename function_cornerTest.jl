function corner_g(y, pred, c_L, c_U, d, x_now, edge, origin, destination)
#     Len = length(d)
#     c_corner = zeros(length(edge[:,1]))
    c_corner = y.*c_U + (1 .- y).*c_L + x_now.*d
#     for i = 1:length(edge[:,1])
#         if y[i] == 1
#             c_corner[i] = c_U[i] + x_now[i]*d[i]
#         else
#             c_corner[i] = c_L[i] + x_now[i]*d[i]
#         end
#     end
#     println("c_corner = ", c_corner)
#     println(edge)
#     println(x_now)
#     println(c_U)
#     println(d)
#     mult = 
    π = zeros(length(edge[:,1])) .+ sum(c_U[k] + x_now[k]*d[k] for k = 1:length(edge[:,1])) .+ 1
#     println(π)
    u = origin
    y_index = findall(y.==1)
    π[u] = 0
    S = [u]
#     println("y_index = ", y_index)
    while u != destination
#         println("u = ", u)
        found_node = false
        count = 0
        for i in y_index
#         while found_node == false 
#             count = count + 1
#             i = y_index[count]
#             println("i = ", i)
            if edge[i,1] == u
#                 found_node = true
                v = edge[i,2]
                π[v] = π[u] + c_corner[i]
                u = v
                push!(S, u)
                break
            end
        end
    end
#     println("S = ", S)
#     println(π)
    B = copy(S) #+ 0
    
    opt = true
    while opt == true && length(B) > 0
#         println("\nB = ", B)
        u = B[1]
        B[1], B[length(B)] = B[length(B)], B[1]
        pop!(B)
        myArcs = findall(edge[:,1].==u)
#         println("empty? ", isempty(myArcs))
        if isempty(myArcs) == false #length(myArcs) > 0
#             println(myArcs)
            for i in myArcs
                if y[i] != 1
                    v = edge[i,2]
                    if π[v] > π[u] + c_corner[i]
                        π[v] = π[u] + c_corner[i]
                        if (v in B) == false #length(find(B.==v)) == 0
                            push!(B, v)
                        end
                        if (v in S) == true #length(find(S.==v)) > 0
                            opt = false
                        end
                    end
                end
            end
        else
            succ = findall(pred.==u)
#             println("succ = ", succ)
            for node in succ
                if (node in B) == true
                    push!(B, node)
                end
            end
#             for i = 1:length(succ)
#                 if length(findall(B.==succ[i])) > 0
#                     push!(B, succ[i])
#                 end
#             end
        end
#         println("opt = ", opt)
    end
    return opt, c_corner
end