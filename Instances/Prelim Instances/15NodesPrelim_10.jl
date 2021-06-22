edge = [1 4; 1 5; 1 7; 2 15; 3 4; 4 5; 4 6; 4 8; 4 9; 4 11; 4 15; 5 7; 6 13; 6 14; 7 10; 8 6; 8 10; 8 14; 8 15; 9 5; 9 8; 9 13; 10 5; 10 15; 11 4; 11 5; 11 12; 11 13; 12 3; 12 8; 12 11; 12 13; 12 15; 13 2; 13 3; 13 7; 13 8; 14 5; 14 8]
cL_orig = [26.0, 40.0, 60.0, 127.0, 8.0, 15.0, 19.0, 41.0, 50.0, 65.0, 108.0, 13.0, 71.0, 79.0, 31.0, 13.0, 23.0, 58.0, 68.0, 36.0, 8.0, 45.0, 55.0, 40.0, 71.0, 55.0, 9.0, 15.0, 94.0, 42.0, 9.0, 12.0, 23.0, 100.0, 104.0, 53.0, 53.0, 87.0, 57.0]
cU_orig = [26.0, 40.0, 60.0, 127.0, 8.0, 15.0, 19.0, 41.0, 50.0, 65.0, 108.0, 31.0, 71.0, 79.0, 31.0, 27.0, 23.0, 68.0, 78.0, 36.0, 8.0, 45.0, 55.0, 54.0, 71.0, 55.0, 9.0, 15.0, 94.0, 42.0, 9.0, 12.0, 37.0, 116.0, 104.0, 71.0, 53.0, 87.0, 57.0]
d = [15.0, 20.0, 16.0, 7.0, 13.0, 13.0, 2.0, 5.0, 17.0, 4.0, 13.0, 10.0, 12.0, 8.0, 17.0, 20.0, 3.0, 5.0, 11.0, 4.0, 17.0, 20.0, 14.0, 17.0, 16.0, 20.0, 14.0, 6.0, 5.0, 18.0, 14.0, 9.0, 18.0, 15.0, 8.0, 20.0, 14.0, 14.0, 20.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
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
delta2 = 3.5
last_node = maximum(edge)
b=2
Grp = "15NodesPrelim"