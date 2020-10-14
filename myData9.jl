
edge = [1 6; 1 7; 1 9; 2 3; 2 4; 2 5; 2 7; 2 9; 3 9; 4 2; 4 3; 4 5; 4 7; 5 1; 5 3; 5 4; 7 1; 7 2; 7 5; 8 1; 8 3; 8 4; 9 5; 9 6]
cL_orig = [31.0, 8.0, 33.0, 1.0, 28.0, 17.0, 1.0, 9.0, 5.0, 9.0, 10.0, 10.0, 20.0, 39.0, 19.0, 32.0, 28.0, 7.0, 20.0, 3.0, 13.0, 12.0, 24.0, 6.0]
cU_orig = [31.0, 38.0, 33.0, 1.0, 28.0, 17.0, 16.0, 9.0, 5.0, 9.0, 10.0, 20.0, 20.0, 39.0, 19.0, 32.0, 28.0, 33.0, 20.0, 3.0, 22.0, 12.0, 24.0, 10.0]
d = [31.0, 27.0, 4.0, 9.0, 40.0, 1.0, 1.0, 7.0, 1.0, 9.0, 9.0, 37.0, 24.0, 14.0, 11.0, 24.0, 1.0, 9.0, 5.0, 8.0, 33.0, 12.0, 14.0, 16.0]


Len = length(cL_orig)

c_orig = 0.5*(cL_orig+cU_orig)

#STARTING SOLUTION:
yy = zeros(Len)

yy[21] = 1
yy[24] = 1
yy[9] = 1

SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)

p = [1.0]
g = [SP_init]
h = [0.0]

origin = 8
destination = 6

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
