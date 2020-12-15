
edge = [1 2; 1 3; 2 3]
cL_orig =  [5.0, 9.0, 5.0]
cU_orig =  [5.0, 13, 5.0]
d =  [1.0, 1.0, 1.0]
Len = length(cL_orig)

c_orig = 0.5*(cL_orig+cU_orig)

#STARTING SOLUTION:
yy = zeros(Len)

yy[2] = 1


SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)

p = [1.0]
g = [SP_init]
h = [0.0]

origin = 1
destination = 3

last_node = maximum(edge)
all_nodes = collect(1:last_node)

M_orig = zeros(Len)

for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end

case = 0
delta1 = 1e-6
delta2 = 0.01
α = 1.0

last_node = maximum(edge)
