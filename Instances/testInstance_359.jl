edge = [1 3; 1 4; 1 9; 1 13; 2 5; 2 8; 2 13; 2 14; 3 6; 3 8; 3 10; 3 11; 3 13; 4 6; 4 12; 4 13; 5 4; 5 9; 5 10; 5 11; 5 13; 5 14; 6 2; 6 4; 6 14; 7 2; 7 9; 7 11; 7 15; 8 2; 8 3; 8 4; 8 7; 9 3; 9 8; 9 11; 9 12; 10 3; 10 4; 10 5; 10 6; 10 8; 10 14; 10 15; 11 2; 11 8; 11 12; 11 13; 11 14; 12 9; 12 11; 12 15; 13 2; 13 6; 13 9; 13 15; 14 2; 14 8; 14 9; 14 11; 14 15]
cL_orig = [33.0, 24.0, 287.0, 371.0, 98.0, 290.0, 213.0, 145.0, 89.0, 59.0, 258.0, 140.0, 419.0, 18.0, 114.0, 346.0, 42.0, 77.0, 135.0, 228.0, 250.0, 40.0, 195.0, 45.0, 198.0, 34.0, 22.0, 30.0, 286.0, 8.0, 136.0, 22.0, 42.0, 110.0, 27.0, 29.0, 141.0, 73.0, 36.0, 84.0, 35.0, 1.0, 170.0, 225.0, 278.0, 69.0, 7.0, 8.0, 107.0, 40.0, 49.0, 13.0, 145.0, 286.0, 125.0, 68.0, 64.0, 252.0, 18.0, 1.0, 9.0]
cU_orig = [49.0, 24.0, 287.0, 371.0, 98.0, 290.0, 213.0, 145.0, 89.0, 263.0, 258.0, 140.0, 419.0, 18.0, 114.0, 346.0, 42.0, 77.0, 135.0, 228.0, 250.0, 40.0, 195.0, 77.0, 288.0, 316.0, 22.0, 30.0, 286.0, 8.0, 136.0, 22.0, 42.0, 110.0, 63.0, 29.0, 141.0, 73.0, 36.0, 84.0, 35.0, 199.0, 170.0, 225.0, 278.0, 69.0, 7.0, 8.0, 155.0, 40.0, 49.0, 21.0, 359.0, 286.0, 125.0, 68.0, 64.0, 252.0, 18.0, 9.0, 9.0]
d = [37.0, 5.0, 29.0, 8.0, 24.0, 39.0, 34.0, 7.0, 32.0, 31.0, 1.0, 23.0, 35.0, 24.0, 2.0, 1.0, 7.0, 9.0, 33.0, 48.0, 39.0, 25.0, 29.0, 28.0, 4.0, 39.0, 47.0, 44.0, 8.0, 21.0, 48.0, 10.0, 41.0, 31.0, 29.0, 30.0, 5.0, 24.0, 49.0, 19.0, 9.0, 38.0, 38.0, 20.0, 40.0, 11.0, 45.0, 29.0, 29.0, 48.0, 46.0, 22.0, 34.0, 10.0, 46.0, 41.0, 24.0, 8.0, 6.0, 30.0, 32.0]
Len = length(d)

yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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