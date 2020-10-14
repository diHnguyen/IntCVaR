
edge = [1 7; 2 1; 3 1; 3 2; 3 5; 3 6; 3 8; 4 1; 4 2; 4 5; 6 2; 6 4; 6 7; 7 6; 7 8; 8 2; 8 6; 8 7]
cL_orig = [61.0, 62.0, 51.0, 21.0, 32.0, 21.0, 90.0, 42.0, 18.0, 7.0, 33.0, 22.0, 13.0, 41.0, 15.0, 29.0, 5.0, 71.0]
cU_orig = [61.0, 67.0, 51.0, 24.0, 84.0, 21.0, 90.0, 42.0, 18.0, 7.0, 34.0, 39.0, 13.0, 77.0, 24.0, 29.0, 8.0, 71.0]
d = [15.0, 38.0, 42.0, 16.0, 1.0, 19.0, 25.0, 2.0, 18.0, 6.0, 13.0, 11.0, 3.0, 12.0, 8.0, 6.0, 3.0, 13.0]

Len = length(cL_orig)

c_orig = 0.5*(cL_orig+cU_orig)

#STARTING SOLUTION:
yy = zeros(Len)

yy[5] = 1

SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)

p = [1.0]
g = [SP_init]
h = [0.0]

origin = 3
destination = 5

last_node = maximum(edge)
all_nodes = collect(1:last_node)

M_orig = zeros(Len)

for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end

case = 0
delta1 = 1e-6
delta2 = 2.0
α = 1.0
last_node = maximum(edge)
