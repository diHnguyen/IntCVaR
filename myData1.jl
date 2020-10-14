
edge = [2 4; 2 5; 3 4; 4 3; 4 5; 5 1]
cL_orig = [58.0, 11.0, 15.0, 45.0, 13.0, 94.0]
cU_orig = [60.0, 95.0, 15.0, 59.0, 26.0, 94.0]
d = [32.0, 31.0, 32.0, 47.0, 44.0, 38.0]

Len = length(cL_orig)

α = 1.0

c_orig = 0.5*(cL_orig+cU_orig)

#STARTING SOLUTION:
yy = zeros(Len)

yy[6] = 1
yy[2] = 1

SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)

p = [1.0]
g = [SP_init]
h = [0.0]

origin = 2
destination = 1

last_node = maximum(edge)
all_nodes = collect(1:last_node)

M_orig = zeros(Len)

for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end

case = 0
delta1 = 1e-6
delta2 = 1.0

last_node = maximum(edge)
