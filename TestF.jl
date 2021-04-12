#Ready for upload
using JuMP 
using TimerOutputs

M = 10
K1 = collect(1:M)
K2 = collect(1:M)
K3 = collect(1:M)

# deleteF = time()

# for k = 1:M
#     deleteat!(K1, K1.==k)
# end
# deleteF = time() - deleteF
# println(deleteF)

filterF = time()
for k = 1:M
    filter!(x->x != k,K2)
    println(K2)
end
filterF = time() - filterF
println(filterF)

# filterF = time()
# for k = 1:M
#     filter!(x->x != k,K2)
# end
# filterF = time() - filterF
# println(filterF)

# deleteF = time()

# for k = 1:M
#     deleteat!(K3, K3.==k)
# end
# deleteF = time() - deleteF

# println(deleteF)