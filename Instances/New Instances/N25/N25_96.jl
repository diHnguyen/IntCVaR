edge = [1 16; 2 25; 3 4; 3 5; 3 10; 3 17; 3 22; 4 11; 4 17; 5 12; 5 15; 5 22; 6 7; 6 19; 7 8; 7 17; 7 24; 8 7; 8 18; 9 15; 10 13; 10 14; 10 23; 10 25; 11 7; 11 13; 11 15; 12 9; 12 14; 13 2; 13 7; 13 12; 14 7; 14 8; 14 10; 14 20; 15 6; 15 21; 16 6; 16 13; 17 8; 17 16; 17 22; 18 10; 18 17; 19 9; 20 6; 20 12; 21 10; 21 11; 21 12; 21 15; 21 18; 22 3; 22 10; 22 11; 23 14; 23 17; 23 18; 23 22; 23 24; 24 3; 24 9; 24 11]
cL_orig = [138.0, 228.0, 10.0, 25.0, 65.0, 135.0, 193.0, 74.0, 126.0, 65.0, 95.0, 172.0, 11.0, 132.0, 0.0, 96.0, 172.0, 10.0, 96.0, 57.0, 21.0, 45.0, 129.0, 147.0, 43.0, 18.0, 39.0, 35.0, 7.0, 101.0, 60.0, 8.0, 71.0, 62.0, 43.0, 61.0, 95.0, 61.0, 98.0, 35.0, 89.0, 15.0, 45.0, 84.0, 6.0, 98.0, 131.0, 81.0, 111.0, 90.0, 86.0, 58.0, 31.0, 184.0, 123.0, 109.0, 80.0, 49.0, 46.0, 4.0, 15.0, 208.0, 148.0, 134.0]
cU_orig = [158.0, 228.0, 10.0, 25.0, 65.0, 135.0, 193.0, 74.0, 126.0, 65.0, 95.0, 172.0, 11.0, 132.0, 20.0, 96.0, 172.0, 10.0, 110.0, 57.0, 35.0, 45.0, 141.0, 147.0, 43.0, 18.0, 39.0, 35.0, 23.0, 119.0, 60.0, 8.0, 71.0, 62.0, 43.0, 61.0, 95.0, 61.0, 98.0, 35.0, 89.0, 15.0, 45.0, 84.0, 6.0, 98.0, 145.0, 81.0, 111.0, 108.0, 86.0, 58.0, 31.0, 204.0, 123.0, 109.0, 90.0, 65.0, 46.0, 20.0, 15.0, 208.0, 148.0, 134.0]
d = [20.0, 4.0, 19.0, 2.0, 12.0, 2.0, 8.0, 5.0, 7.0, 3.0, 15.0, 11.0, 3.0, 19.0, 15.0, 5.0, 13.0, 10.0, 15.0, 16.0, 19.0, 7.0, 7.0, 14.0, 1.0, 9.0, 8.0, 6.0, 16.0, 9.0, 6.0, 2.0, 18.0, 20.0, 6.0, 6.0, 6.0, 5.0, 18.0, 12.0, 1.0, 11.0, 8.0, 14.0, 14.0, 7.0, 5.0, 2.0, 6.0, 18.0, 6.0, 18.0, 17.0, 11.0, 16.0, 13.0, 16.0, 14.0, 18.0, 1.0, 20.0, 14.0, 3.0, 19.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
β = 0.166