edge = [1 3; 1 5; 2 8; 2 10; 3 4; 3 5; 3 6; 3 7; 3 9; 3 11; 4 3; 4 9; 5 2; 5 13; 6 5; 6 13; 7 5; 7 12; 8 13; 8 14; 9 3; 9 4; 9 11; 10 3; 10 6; 10 8; 10 9; 10 14; 10 15; 11 5; 11 6; 11 12; 11 13; 12 6; 12 10; 13 4; 14 3; 14 13]
cL_orig = [21.0, 36.0, 64.0, 81.0, 8.0, 15.0, 23.0, 39.0, 62.0, 85.0, 0.0, 50.0, 29.0, 82.0, 12.0, 70.0, 10.0, 54.0, 48.0, 57.0, 57.0, 52.0, 18.0, 69.0, 41.0, 22.0, 11.0, 45.0, 55.0, 65.0, 44.0, 8.0, 18.0, 54.0, 12.0, 94.0, 110.0, 12.0]
cU_orig = [21.0, 36.0, 64.0, 81.0, 20.0, 31.0, 35.0, 39.0, 62.0, 85.0, 20.0, 50.0, 29.0, 82.0, 12.0, 70.0, 22.0, 54.0, 48.0, 57.0, 57.0, 52.0, 18.0, 69.0, 41.0, 22.0, 11.0, 45.0, 55.0, 65.0, 54.0, 8.0, 18.0, 74.0, 28.0, 94.0, 110.0, 12.0]
d = [15.0, 18.0, 5.0, 12.0, 11.0, 2.0, 11.0, 16.0, 17.0, 19.0, 15.0, 20.0, 7.0, 8.0, 2.0, 15.0, 6.0, 8.0, 13.0, 17.0, 15.0, 16.0, 17.0, 19.0, 9.0, 3.0, 9.0, 3.0, 9.0, 14.0, 8.0, 1.0, 18.0, 4.0, 11.0, 14.0, 20.0, 12.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0]
SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)
p = [1.0]
g = [SP_init]
h = [0.0]
origin = 1
destination = 15
last_node = maximum(edge)
all_nodes = collect(1:last_node)
M_orig = zeros(Len)
for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end
delta1 = 1.0
delta2 = 2.0
last_node = maximum(edge)
b=2
Grp = "15NodesPrelim"