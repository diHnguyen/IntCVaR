edge = [1 2; 1 12; 1 20; 2 5; 3 5; 3 12; 4 11; 4 18; 5 6; 6 19; 7 14; 8 18; 9 5; 10 11; 10 13; 10 14; 11 18; 12 2; 12 6; 12 17; 13 8; 13 10; 13 20; 14 7; 15 4; 15 8; 15 13; 15 14; 15 19; 16 3; 16 14; 17 3; 17 4; 17 15; 17 19; 18 13; 18 16; 19 9; 19 15]
cL_orig = Any[15.0, 106.0, 190.0, 31.0, 10.0, 92.0, 65.0, 131.0, 7.0, 123.0, 68.0, 99.0, 45.0, 0.0, 25.0, 30.0, 70.0, 99.0, 55.0, 42.0, 50.0, 29.0, 72.0, 73.0, 98.0, 66.0, 8.0, 2.0, 31.0, 119.0, 14.0, 138.0, 126.0, 15.0, 12.0, 47.0, 12.0, 104.0, 38.0]
cU_orig = Any[15.0, 106.0, 190.0, 31.0, 30.0, 92.0, 65.0, 147.0, 7.0, 135.0, 68.0, 99.0, 45.0, 16.0, 25.0, 46.0, 70.0, 99.0, 55.0, 60.0, 50.0, 29.0, 72.0, 73.0, 114.0, 66.0, 26.0, 20.0, 41.0, 135.0, 26.0, 138.0, 126.0, 33.0, 30.0, 57.0, 24.0, 104.0, 38.0]
d = Any[1.0, 4.0, 10.0, 5.0, 9.0, 1.0, 11.0, 9.0, 18.0, 15.0, 9.0, 16.0, 9.0, 8.0, 1.0, 3.0, 19.0, 14.0, 19.0, 7.0, 9.0, 5.0, 4.0, 6.0, 15.0, 18.0, 2.0, 6.0, 14.0, 20.0, 19.0, 8.0, 3.0, 16.0, 8.0, 2.0, 5.0, 9.0, 12.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
delta2 = 3.0
b = 2
last_node = maximum(edge)
Grp = "0.1Density"