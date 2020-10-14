
edge = [1 4; 1 5; 1 6; 3 1; 3 4; 3 6; 4 1; 4 5; 4 6; 5 4; 5 6; 6 1; 6 2]
cL_orig = [16.0, 83.0, 35.0, 84.0, 55.0, 39.0, 87.0, 23.0, 34.0, 21.0, 5.0, 32.0, 93.0]
cU_orig = [40.0, 83.0, 35.0, 84.0, 55.0, 73.0, 87.0, 41.0, 34.0, 72.0, 5.0, 76.0, 93.0]
d = [17.0, 75.0, 20.0, 55.0, 49.0, 20.0, 52.0, 20.0, 32.0, 36.0, 3.0, 71.0, 14.0]

Len = length(cL_orig)

c_orig = 0.5*(cL_orig+cU_orig)

#STARTING SOLUTION:
yy = zeros(Len)

yy[6] = 1
yy[13] = 1

SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)

p = [1.0]
g = [SP_init]
h = [0.0]

origin = 3
destination = 2

last_node = maximum(edge)
all_nodes = collect(1:last_node)

M_orig = zeros(Len)

for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end

case = 0
delta1 = 1e-6
delta2 = 3.5
α = 1.0
last_node = maximum(edge)
