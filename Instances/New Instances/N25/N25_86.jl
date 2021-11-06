edge = [1 3; 2 3; 2 5; 2 7; 2 12; 3 10; 4 16; 4 19; 5 2; 5 12; 6 20; 7 22; 8 14; 8 19; 9 13; 9 15; 10 14; 10 19; 10 20; 10 21; 10 22; 11 9; 11 18; 11 20; 11 23; 12 6; 12 11; 12 24; 13 11; 14 4; 14 8; 14 21; 14 23; 15 5; 15 8; 15 18; 16 8; 16 17; 16 19; 16 22; 16 24; 17 2; 17 10; 17 16; 17 19; 17 21; 18 7; 18 9; 18 19; 19 12; 19 25; 20 19; 21 6; 21 12; 21 19; 22 10; 22 23; 23 3; 23 4; 23 14; 23 25; 24 21]
cL_orig = [24.0, 11.0, 26.0, 55.0, 105.0, 74.0, 120.0, 151.0, 35.0, 75.0, 145.0, 147.0, 58.0, 107.0, 37.0, 57.0, 30.0, 77.0, 105.0, 105.0, 124.0, 10.0, 69.0, 90.0, 119.0, 47.0, 9.0, 107.0, 22.0, 91.0, 61.0, 62.0, 95.0, 88.0, 71.0, 27.0, 83.0, 10.0, 29.0, 60.0, 75.0, 149.0, 70.0, 12.0, 16.0, 37.0, 105.0, 91.0, 9.0, 73.0, 57.0, 7.0, 147.0, 88.0, 13.0, 123.0, 9.0, 201.0, 190.0, 87.0, 16.0, 32.0]
cU_orig = [24.0, 11.0, 26.0, 55.0, 105.0, 74.0, 120.0, 151.0, 35.0, 75.0, 145.0, 147.0, 58.0, 107.0, 37.0, 57.0, 42.0, 93.0, 105.0, 123.0, 124.0, 28.0, 69.0, 90.0, 119.0, 67.0, 9.0, 123.0, 22.0, 111.0, 61.0, 80.0, 95.0, 106.0, 71.0, 27.0, 83.0, 10.0, 29.0, 60.0, 75.0, 149.0, 70.0, 12.0, 34.0, 37.0, 119.0, 91.0, 9.0, 73.0, 57.0, 7.0, 147.0, 88.0, 27.0, 123.0, 9.0, 201.0, 190.0, 87.0, 16.0, 32.0]
d = [14.0, 7.0, 6.0, 16.0, 8.0, 6.0, 2.0, 9.0, 6.0, 7.0, 14.0, 9.0, 15.0, 10.0, 20.0, 1.0, 9.0, 5.0, 15.0, 13.0, 5.0, 8.0, 2.0, 18.0, 3.0, 20.0, 2.0, 19.0, 4.0, 10.0, 3.0, 4.0, 10.0, 2.0, 1.0, 14.0, 2.0, 10.0, 3.0, 14.0, 13.0, 8.0, 15.0, 2.0, 3.0, 15.0, 20.0, 6.0, 14.0, 1.0, 4.0, 6.0, 19.0, 11.0, 15.0, 1.0, 3.0, 13.0, 8.0, 12.0, 13.0, 5.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
β = 0.78