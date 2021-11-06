edge = [1 3; 1 4; 2 7; 2 9; 2 12; 2 15; 2 22; 2 25; 3 23; 4 6; 5 7; 5 9; 5 15; 5 18; 6 10; 6 11; 7 3; 7 21; 8 12; 8 13; 8 19; 8 20; 9 6; 10 3; 10 9; 10 21; 11 12; 11 24; 12 4; 12 16; 13 14; 13 18; 14 22; 15 7; 15 23; 16 20; 17 5; 17 10; 17 15; 18 2; 18 12; 18 17; 19 10; 19 13; 19 18; 20 22; 20 24; 21 11; 21 12; 22 8; 22 16; 22 17; 23 8; 23 18; 24 8; 24 15; 24 17]
cL_orig = [17.0, 31.0, 52.0, 70.0, 105.0, 131.0, 197.0, 226.0, 204.0, 15.0, 17.0, 38.0, 105.0, 128.0, 42.0, 55.0, 34.0, 138.0, 45.0, 49.0, 107.0, 110.0, 15.0, 72.0, 11.0, 114.0, 5.0, 126.0, 82.0, 29.0, 15.0, 47.0, 78.0, 76.0, 67.0, 35.0, 118.0, 69.0, 18.0, 162.0, 46.0, 13.0, 95.0, 60.0, 9.0, 22.0, 36.0, 94.0, 87.0, 144.0, 60.0, 48.0, 148.0, 49.0, 158.0, 92.0, 75.0]
cU_orig = [31.0, 31.0, 52.0, 70.0, 105.0, 131.0, 209.0, 226.0, 204.0, 15.0, 17.0, 38.0, 105.0, 128.0, 42.0, 55.0, 52.0, 138.0, 45.0, 49.0, 107.0, 126.0, 35.0, 72.0, 11.0, 114.0, 5.0, 126.0, 82.0, 49.0, 15.0, 47.0, 78.0, 76.0, 85.0, 35.0, 118.0, 81.0, 18.0, 162.0, 64.0, 13.0, 95.0, 60.0, 9.0, 22.0, 36.0, 112.0, 103.0, 144.0, 60.0, 48.0, 148.0, 49.0, 158.0, 92.0, 75.0]
d = [20.0, 2.0, 13.0, 18.0, 17.0, 6.0, 4.0, 13.0, 18.0, 3.0, 13.0, 12.0, 11.0, 1.0, 18.0, 10.0, 8.0, 17.0, 17.0, 13.0, 3.0, 8.0, 16.0, 3.0, 11.0, 13.0, 12.0, 19.0, 8.0, 12.0, 16.0, 18.0, 11.0, 3.0, 10.0, 11.0, 9.0, 9.0, 4.0, 12.0, 11.0, 4.0, 15.0, 16.0, 16.0, 15.0, 10.0, 9.0, 18.0, 20.0, 10.0, 16.0, 13.0, 5.0, 15.0, 20.0, 5.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)
p = [1.0]
g = [SP_init]
h = [0.0]
origin = 1
destination = 25
last_node = maximum(edge)
all_nodes = collect(1:last_node)
M_orig = zeros(Len)
for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end
delta1 = 1.0
delta2 = 2.0
b = 7
last_node = maximum(edge)
β = 0.379
