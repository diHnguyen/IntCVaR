using DataFrames
using Polynomials
using JuMP
# using TimerOutputs
# using Plots

function findRange((l1,u1),(l2,u2)) #Used in Convolution
    r = u1-l1
    R = u2-l2
    bounds = Array{Float64}(undef, (0, 2))
    if R != r
        bounds = [0 r; r R; R R+r]
#         bounds =[ ]
        bounds = bounds.+(l1+l2)
    else
        bounds = [0 r; r r+R]
        bounds = bounds.+(l1+l2)
    end
    return bounds
end
function timeShift(poly, direction, r) #Used in Convolution
    c = coeffs(poly)
    myRoot = 0
    if direction == "left"
        myRoot = -r
    end
    
    if direction == "right"
        myRoot = r
    end
    multiplier_poly = Polynomial([1]) 
    # println(multiplier_poly)
    my_poly = Polynomial([0])

    for i = 1:length(c)
        my_poly = my_poly + multiplier_poly*c[i]
        multiplier_poly = multiplier_poly*(fromroots([myRoot])) 
    end
    return my_poly
end

function Convolve(cL, cU, Len, unc_arcs)
#     cL = [0,0,0].+1
#     cU = [3,1,2].+1
#     println("cL = ", cL)
#     println("cU = ", cU)
#     global Len
    #Temporary Creating ranges here:
    df_Uniform = DataFrame(LB = Float64[], UB = Float64[], r = Float64[])
    
    
#     println("Deterministic Arcs = ", det_arcs)
    for i = 2:length(unc_arcs)
        k = unc_arcs[i]
        push!(df_Uniform,(cL[k],cU[k],cU[k]-cL[k]))
    end
#     println("dfUniform ", df_Uniform)
    numOrder = length(cL)  #maxOrder = numOrder - 1 since we need to account for degree 0
    df = DataFrame(PolyNum = String[], Poly = Polynomial{Float64}[], l = Float64[], u = Float64[], leftShift = Float64[])
    first_arc = unc_arcs[1]
    push!(df,("", 1/(cU[first_arc]-cL[first_arc]), cL[first_arc],cU[first_arc], 0.0))
    numPoly = 1
    ii=0
    
    while length(df_Uniform[:,:LB]) > 0
#         println(df)

        l1 = df_Uniform[1,:LB]
        u1 = df_Uniform[1,:UB]
        r1 = df_Uniform[1,:r]
        
        poly1 = Polynomial([1/r1])
        
        #Time shift Uniform Var.###
        l1LeftShift = l1
        ####Update (l1,u1)
        u1 = r1
        l1 = 0
        
        poly_to_Integrate = length(df[:, :PolyNum])
        ii=ii+1
       
        for k = 1:poly_to_Integrate #numPoly #numOrder-1
            l2 = df[k,:l]
            u2 = df[k,:u]
            r2 = u2 - l2
            
            poly2 = df[k,:Poly]
            
            #Shift poly2:
            poly2 = timeShift(poly2, "left", l2)
            totalLeftShift = l1LeftShift + l2 + df[k,:leftShift]
            u2 = r2
            l2 = 0
            
            bounds = []
            int_L = 0
            int_U = 0
            
            abs_L = l2
            abs_U = u2

            if r2 > r1
                bounds = findRange((l1,u1), (l2,u2))
            else
                bounds = findRange((l2,u2), (l1,u1))
            end

            #This returns the function that we need to integrate wrt different sets of bounds
            p_combine = poly1 * poly2
            
            #Integrate p_combine: Same for all 2 or 3 ranges below
            p_int = integrate(p_combine)

            #Calculate the lower-bound polynomial with variable 
            #E.g. if the integration goes from t-4 to t then the lower bound is evaluated at t-4
            #E.g. if the integration goes from t-3 to t-1 then the upper bound is evaluated at t-1
            #Note that for the lower bound integral we use fromroots([int_U]) since t-int_U gives us a lower bound
            #Note that for the upper bound integral we use fromroots([int_L]) since t-int_L gives us an upper bound
            
            #Evaluate p_int at the lower bound t-root_L is
            #similar to perform time shift to the right or when = "post" 
            
            my_poly_L = timeShift(p_int, "right", r1)
            
            #Evaluate with 1st resulting range: from L to t
            L = bounds[1,1]
            U = bounds[1,2]
            p_1 = p_int - p_int(l2) #Note that p_combine evaluated at t returns the same polynomial 
            df[k,:] = df[k,:PolyNum]*".1", p_1, L, U, totalLeftShift         

            #Deal with 2nd resulting range: from t-some_range to t
            #If there are only 2 resulting ranges then we don't need to calculate this.
            row = 2
            if length(bounds[:,1]) > 2
                numPoly = numPoly + 1
                p_2 = Polynomial([0])
                if r2 > r1
                    p_2 = p_int - my_poly_L
                end
            
                if r2 < r1
                    interval_2 = p_int(u2) - p_int(l2)
                    p_2 = Polynomial([interval_2])
                end
                push!(df,(df[k,:PolyNum]*"."*string(row), p_2, bounds[row,1], bounds[row,2],totalLeftShift))
                row = row + 1
            end
            
            p_3 = p_int(u2) - my_poly_L
            numPoly = numPoly + 1
            push!(df,(df[k,:PolyNum]*"."*string(row),p_3,bounds[row,1],bounds[row,2],totalLeftShift))

        end
        delete!(df_Uniform, 1)
        
    end
    return df, numPoly
end

function FindCVaR(β,α,gL,gU,df_cellPoly, pCell)
    tol = 1e-3
    W = -1
    V = 0
    α_L = gL
    α_U = gU
    #The original
    VaR = α
#     println("1-β = ", 1-β)
#     println("1-β - W = ", 1-β-W)
    iter = 1
    lastVaR = -0.5
    W_k = zeros(3)
    while abs(1-β - W) > tol
        
        println("\nVaR Guess = ", VaR)
        println("W current = ", W)
        W = 0
        
        for k=1:3 #k in K
#             println()
#             println("lastVaR = ", lastVaR)
#             cL = [0,0]
#             cU = [4,2]
#             println("k = ", k)
#             df, numPoly = Convolve(cL, cU)
            df = df_cellPoly[k,:df]
#             println(df)
            numPoly = df_cellPoly[k,:NUMPOLY]
            det_Shift = df_cellPoly[k,:DETSHIFT]
#             df[:,:w] = zeros(numPoly) #.*(-1)
#             println(df)
#             break
            W_k[k] = 0
#             println("df_w  = ", df[:, :w])
            for i = 1:numPoly
                r_l = df[i,:l]
                r_u = df[i,:u]
                rightShift = df[i,:leftShift]
                
                poly = df[i,:Poly]
                p_poly = integrate(poly)
#                 println("entire Prob = ", p_poly(r_u) - p_poly(r_l))
                
#                 println("r_l = ", r_l)
#                 println("r_u = ", r_u)
#                 println(i,". ", df[i,:])
                w_i = 0
                if r_l + rightShift + det_Shift <= VaR  #df[i,:leftShift]
#                     println()
#                     println("r_l + rightShift + det_Shift = ", r_l + rightShift + det_Shift) 
#                     println("r_u + rightShift + det_Shift = ", r_u + rightShift + det_Shift)
    #                 if df[i,:w] < 0 || r_u > VaR - rightShift
#                     println("lastVaR = ", lastVaR)
                    if (lastVaR < r_u + rightShift + det_Shift) || (VaR < r_u + rightShift+ det_Shift)

#                         println("r_u = ", r_u)
#                         println("VaR - rightShift - det_Shift = ", VaR-rightShift- det_Shift)
                        poly = df[i,:Poly]
                        p_poly = integrate(poly)
                        
                        u = min(VaR-rightShift-det_Shift, r_u) 
#                         println("Integrating from ", r_l ," to ", u)
    #                     l = df[i,:l]
                        w_i = p_poly(u) - p_poly(r_l)
#                         println(i, ". w_i = ", w_i)
                        df[i,:w] = w_i
                    else
                        w_i = df[i,:w]
                    end
                    
                    
                    
    #                 if r_u > VaR - rightShift
    #                     df[i,:w] = -1
    #                 end
                end
                
                W_k[k] = W_k[k] + w_i
                
            end
#             println(k,". W_k = ", W_k,"\n")
            W = W + W_k[k]*pCell[k]
#             println("Probability = ", W) 
        end
        println("W = ", W)
        lastVaR = VaR
        if VaR > α_L && W <= 1-β
            α_L = VaR
        end
        if VaR < α_U && W >= 1-β
            α_U = VaR
        end
        if abs(1-β - W) > tol
           
            VaR = (α_U + α_L)/2
        end
#         println("1-β - W = ", 1-β-W)
#         iter = iter +1
#         if iter > 5
#             break
#         end
    end
    println("Final β-VaR = ", VaR)
    println("1-β = ", 1-β)
    println("Final W = ", W)
    println("W_k = ", W_k)
    V=0
    total_p = 0
    V_k = zeros(3)
#     VW_k = zeros(3)
    for k = 1:3 #Replace with the size of K here
        if W_k[k] > 0
            df = df_cellPoly[k,:df]
            numPoly = df_cellPoly[k,:NUMPOLY]
            det_Shift = df_cellPoly[k,:DETSHIFT]
            V = 0
            for i = 1:numPoly
                r_l = df[i,:l]
                r_u = df[i,:u]
                rightShift = df[i,:leftShift]
                poly = df[i,:Poly]
                if r_l + rightShift + det_Shift <= VaR 
                    println("Assessing V, poly ", i)
                    u = min(VaR-rightShift-det_Shift, r_u) 
                    println(poly)
                    println(fromroots([-(rightShift+det_Shift)]))
                    
                    
                    e_poly = integrate(poly*fromroots([-(rightShift+det_Shift)]))
                    v_i = e_poly(u) - e_poly(r_l)
    #                 println("v_i = ", v_i)
                    V = V + v_i #*df[i,:w]
                end
        
            end
            V_k[k] = V
            println(k,". V_k = ", V_k[k], "; W_k = ", W_k[k], "; p_k = ", pCell[k])
#             println()
#             V = V + V_k/W_k[k]*pCell[k]
#             total_p = total_p + pCell[k]
        end
    end
#     println("Final V = ", V)
    println("β-CVaR = ", sum(V_k[k]*pCell[k] for k =1:3)/sum(W_k[k]*pCell[k] for k =1:3)) #/W)
end

# β = 0.5
# α = 1
# gL = 0
# gU = 6
# FindCVaR(1-β,α,gL,gU)


#Verify with Monte Carlo 
# R = 1000000
# for rep = 1:5
#     count = 0
#     avg = 0
#     for i = 1:R
#         a = rand()*4
#         b = rand()*2

#         if a + b <= 3
#             count = count + 1
#             avg = avg + a+b
#         end
#     end
#     println("Monte Carlo: ", avg/count) 
# end