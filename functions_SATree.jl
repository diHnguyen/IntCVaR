include("./functions_sharedTreePath.jl")
function find_basic_Delta(myNode, Delta_Final_Pos, Delta_Final_Neg, sign, basic, i_basic, nonbasic, i_nonbasic)
    myNonTreeArcList = []
    myTreeArcList = []
    Delta_List = []
    delta = -1
#     println("myNode = ", myNode)
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
        
        e = i_nonbasic[myNonTreeArcList[1]]
#         println("Arc ", e)
        delta = Delta_Final_Neg[e]
#         println("New delta = ", delta)
        for i = 2 : length(myNonTreeArcList)
            e = i_nonbasic[myNonTreeArcList[i]]
#             println("Arc ", e)
            if delta > Delta_Final_Neg[e] && Delta_Final_Neg[e] >= 0
                delta = Delta_Final_Neg[e]
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
            e = i_basic[myTreeArcList[1]]
            delta = Delta_List[e]
#             println("New delta = ", delta)
        end
        for i = 1 : length(myTreeArcList)
            e = i_basic[myTreeArcList[i]]
#             println("Arc ", e)
            if delta > Delta_List[e] && Delta_List[e] >= 0
                delta = Delta_List[e]
#                 println("New delta = ", delta)
            end
        end
    end
    return delta 
end
