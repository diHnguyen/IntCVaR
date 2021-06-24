edge = [1 6; 1 16; 1 20; 4 17; 5 4; 5 8; 5 13; 5 15; 6 9; 6 17; 7 18; 8 9; 8 17; 9 18; 12 4; 12 5; 12 6; 12 13; 12 17; 13 4; 13 12; 15 4; 15 5; 15 7; 15 13; 15 17; 15 20; 16 6; 16 7; 16 8; 16 17; 17 15; 17 20; 18 5; 18 6; 18 20; 19 17]
cL_orig = Any[39.0, 155.0, 191.0, 125.0, 9.0, 18.0, 71.0, 102.0, 29.0, 114.0, 102.0, 15.0, 86.0, 88.0, 73.0, 58.0, 54.0, 6.0, 51.0, 87.0, 7.0, 108.0, 95.0, 78.0, 22.0, 16.0, 48.0, 99.0, 87.0, 82.0, 15.0, 25.0, 22.0, 131.0, 115.0, 21.0, 19.0]
cU_orig = Any[53.0, 155.0, 191.0, 141.0, 9.0, 32.0, 89.0, 102.0, 39.0, 114.0, 114.0, 15.0, 98.0, 102.0, 93.0, 76.0, 64.0, 6.0, 51.0, 87.0, 7.0, 108.0, 95.0, 78.0, 22.0, 26.0, 48.0, 99.0, 87.0, 82.0, 15.0, 25.0, 38.0, 131.0, 115.0, 21.0, 19.0]
d = Any[16.0, 1.0, 4.0, 13.0, 16.0, 1.0, 17.0, 17.0, 5.0, 17.0, 3.0, 19.0, 18.0, 12.0, 10.0, 1.0, 1.0, 2.0, 13.0, 15.0, 8.0, 15.0, 3.0, 11.0, 10.0, 6.0, 4.0, 15.0, 17.0, 2.0, 11.0, 7.0, 14.0, 18.0, 11.0, 20.0, 13.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]
SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)
p = [1.0]
g = [SP_init]
h = [0.0]
origin = 1
destination = 20
last_node = maximum(edge)
all_nodes = collect(1:last_node)
M_orig = zeros(Len)
for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end
delta1 = 1.0
delta2 = 5.0
b = 2
last_node = maximum(edge)
Grp = "-5Nodes"