
edge = [1 2; 1 3; 2 6; 3 4; 3 5; 4 6; 5 6; 1 7; 7 6]
cL_orig = [2.0, 1, 0, 0, 0, 12, 0, 1, 2]
cU_orig = [2.0, 1,10, 0, 0, 12,20, 1, 6]
d = [100.0, 100.0, 0, 0, 0, 0, 0,100, 0]

global Len = length(cL_orig)

c_orig = 0.5*(cL_orig+cU_orig)

#STARTING SOLUTION:
yy = zeros(Len)

yy[1] = 1
yy[3] = 1


SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)

p = [1.0]

origin = 1
destination = 6

last_node = maximum(edge)
all_nodes = collect(1:last_node)

M_orig = cU_orig - cL_orig

last_node = maximum(edge)
all_nodes = collect(1:last_node)

delta1 = 1e-6
delta2 = 1.0

β = 0.8
b=2