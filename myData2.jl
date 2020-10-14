
edge = [1 3; 1 10; 1 12; 2 3; 2 16; 3 13; 3 14; 4 1; 4 6; 4 9; 4 16; 5 1; 5 6; 6 1; 6 11; 7 14; 7 17; 8 7; 10 2; 10 8; 10 13; 10 17; 11 2; 11 3; 11 8; 12 13; 13 5; 13 9; 14 4; 14 7; 15 1; 15 4; 15 10; 15 12; 15 14; 16 5; 16 7; 17 6; 17 9]
cL_orig = [5.0, 6.0, 6.0, 1.0, 7.0, 4.0, 4.0, 2.0, 6.0, 5.0, 5.0, 2.0, 5.0, 2.0, 5.0, 3.0, 5.0, 4.0, 1.0, 1.0, 7.0, 2.0, 4.0, 2.0, 1.0, 5.0, 4.0, 4.0, 5.0, 2.0, 6.0, 2.0, 2.0, 4.0, 3.0, 2.0, 4.0, 1.0, 1.0]
cU_orig = [5.0, 6.0, 6.0, 1.0, 7.0, 5.0, 4.0, 3.0, 6.0, 5.0, 5.0, 2.0, 5.0, 6.0, 5.0, 3.0, 5.0, 4.0, 1.0, 1.0, 7.0, 2.0, 4.0, 6.0, 1.0, 5.0, 4.0, 4.0, 5.0, 7.0, 6.0, 2.0, 5.0, 6.0, 3.0, 5.0, 4.0, 1.0, 3.0]
d = [3.0, 6.0, 8.0, 9.0, 1.0, 2.0, 4.0, 2.0, 3.0, 5.0, 2.0, 9.0, 4.0, 6.0, 10.0, 9.0, 5.0, 9.0, 9.0, 4.0, 4.0, 6.0, 2.0, 4.0, 8.0, 1.0, 9.0, 7.0, 10.0, 2.0, 9.0, 7.0, 1.0, 6.0, 7.0, 9.0, 9.0, 4.0, 2.0]

Len = length(cL_orig)

c_orig = 0.5*(cL_orig+cU_orig)

#STARTING SOLUTION:
yy = zeros(Len)

yy[34] = 1
yy[26] = 1
yy[28] = 1

SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)

p = [1.0]
g = [SP_init]
h = [0.0]

origin = 15
destination = 9

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
