
edge = [1 7; 2 3; 2 5; 3 6; 4 1; 5 1; 6 4]
cL_orig =  [14, 85, 41, 22, 75, 53,  46.0]
cU_orig =  [18, 85, 74, 27, 79, 53,  69.0]
d =  [13, 22, 67, 5, 2, 44, 14.0]
Len = length(cL_orig)

c_orig = 0.5*(cL_orig+cU_orig)

#STARTING SOLUTION:
yy = zeros(Len)

yy[3] = 1
yy[6] = 1
yy[1] = 1

SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)

p = [1.0]
g = [SP_init]
h = [0.0]

origin = 2
destination = 7

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
