edge = [1 2; 1 3; 1 4; 1 6; 1 10; 1 11; 1 12; 2 10; 2 14; 3 2; 3 5; 3 8; 3 15; 4 2; 4 13; 5 8; 5 10; 5 12; 6 7; 6 8; 6 10; 6 12; 6 13; 7 13; 8 2; 8 4; 8 5; 8 6; 8 7; 8 11; 8 12; 8 13; 8 14; 9 8; 9 11; 9 13; 10 2; 10 3; 10 4; 10 5; 10 8; 11 3; 11 4; 11 7; 11 12; 11 14; 12 2; 12 5; 12 6; 13 6; 13 7; 13 9; 13 10; 13 12; 13 15; 14 2; 14 3; 14 10; 14 13]
cL_orig = [18.0, 16.0, 142.0, 144.0, 78.0, 339.0, 399.0, 184.0, 6.0, 25.0, 9.0, 87.0, 455.0, 24.0, 280.0, 14.0, 60.0, 66.0, 12.0, 12.0, 60.0, 170.0, 218.0, 200.0, 62.0, 94.0, 68.0, 19.0, 1.0, 1.0, 44.0, 125.0, 1.0, 8.0, 0.0, 131.0, 70.0, 283.0, 121.0, 147.0, 21.0, 110.0, 239.0, 2.0, 27.0, 74.0, 356.0, 312.0, 159.0, 87.0, 199.0, 188.0, 115.0, 41.0, 35.0, 339.0, 11.0, 90.0, 3.0]
cU_orig = [32.0, 16.0, 142.0, 144.0, 78.0, 339.0, 399.0, 548.0, 6.0, 25.0, 9.0, 87.0, 455.0, 34.0, 280.0, 14.0, 60.0, 66.0, 64.0, 38.0, 216.0, 170.0, 218.0, 298.0, 118.0, 194.0, 68.0, 19.0, 69.0, 1.0, 44.0, 125.0, 5.0, 8.0, 20.0, 131.0, 70.0, 283.0, 121.0, 147.0, 21.0, 110.0, 239.0, 8.0, 69.0, 74.0, 356.0, 312.0, 159.0, 87.0, 199.0, 188.0, 115.0, 41.0, 35.0, 339.0, 101.0, 90.0, 45.0]
d = [46.0, 46.0, 35.0, 22.0, 5.0, 46.0, 47.0, 35.0, 1.0, 36.0, 7.0, 45.0, 6.0, 16.0, 50.0, 31.0, 12.0, 22.0, 8.0, 21.0, 13.0, 37.0, 16.0, 42.0, 42.0, 21.0, 23.0, 13.0, 3.0, 2.0, 24.0, 30.0, 29.0, 12.0, 12.0, 35.0, 26.0, 13.0, 24.0, 30.0, 30.0, 12.0, 22.0, 20.0, 46.0, 1.0, 33.0, 3.0, 1.0, 47.0, 28.0, 37.0, 35.0, 49.0, 36.0, 6.0, 21.0, 8.0, 23.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1]

c_orig = 0.5*(cL_orig+cU_orig)

SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)


p = [1.0]

g = [SP_init]

h = [0.0]


origin = 1

destination =15

last_node = maximum(edge)
all_nodes = collect(1:last_node)

M_orig = zeros(Len)

for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end

case = 0
delta1 = 1e-6
delta2 = 5
last_node = maximum(edge)