
edge = [1 2; 2 3; 4 3; 4 5; 5 1]
cL_orig = [43.0, 38.0, 19.0, 47.0, 27.0]
cU_orig = [43.0, 57.0, 98.0, 47.0, 80.0]
d = [2.0, 5.0, 83.0, 9.0, 3.0]
Len = length(cL_orig)

c_orig = 0.5*(cL_orig+cU_orig)

#STARTING SOLUTION:
yy = zeros(Len)

yy[3] = 1

SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)

p = [1.0]
g = [SP_init]
h = [0.0]

origin = 4
destination = 3

last_node = maximum(edge)
all_nodes = collect(1:last_node)

M_orig = zeros(Len)

for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end

case = 0
delta1 = 1e-6
delta2 = 1.0
α = 1.0
last_node = maximum(edge)
