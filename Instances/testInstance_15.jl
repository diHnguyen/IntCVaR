edge = [1 5; 1 6; 1 8; 1 10; 1 11; 1 13; 2 8; 2 10; 2 11; 3 6; 3 9; 3 10; 3 12; 3 15; 4 6; 4 9; 4 11; 4 12; 4 13; 5 4; 6 2; 6 8; 6 10; 6 11; 7 3; 7 5; 7 8; 7 15; 8 4; 8 6; 8 7; 8 13; 9 5; 9 10; 10 3; 10 7; 10 9; 10 15; 11 2; 11 7; 11 8; 11 10; 11 14; 11 15; 12 4; 12 5; 12 6; 13 2; 13 7; 13 15; 14 2; 14 3; 14 5; 14 11; 14 13]
cL_orig = [109.0, 60.0, 70.0, 107.0, 193.0, 116.0, 35.0, 11.0, 288.0, 42.0, 67.0, 321.0, 59.0, 320.0, 49.0, 68.0, 249.0, 258.0, 136.0, 40.0, 78.0, 12.0, 40.0, 95.0, 28.0, 87.0, 26.0, 343.0, 180.0, 89.0, 6.0, 134.0, 101.0, 47.0, 3.0, 51.0, 12.0, 78.0, 325.0, 141.0, 57.0, 33.0, 115.0, 178.0, 210.0, 59.0, 151.0, 503.0, 9.0, 67.0, 418.0, 170.0, 323.0, 105.0, 5.0]
cU_orig = [109.0, 82.0, 70.0, 189.0, 193.0, 116.0, 35.0, 11.0, 504.0, 94.0, 67.0, 321.0, 59.0, 324.0, 143.0, 68.0, 249.0, 258.0, 136.0, 40.0, 78.0, 12.0, 326.0, 95.0, 122.0, 87.0, 26.0, 343.0, 180.0, 89.0, 6.0, 134.0, 101.0, 47.0, 3.0, 51.0, 12.0, 242.0, 325.0, 225.0, 87.0, 33.0, 115.0, 178.0, 210.0, 59.0, 151.0, 503.0, 9.0, 67.0, 418.0, 170.0, 323.0, 129.0, 5.0]
d = [15.0, 41.0, 12.0, 3.0, 23.0, 44.0, 39.0, 9.0, 47.0, 25.0, 5.0, 43.0, 39.0, 22.0, 39.0, 50.0, 29.0, 20.0, 13.0, 15.0, 11.0, 33.0, 12.0, 15.0, 47.0, 1.0, 20.0, 3.0, 2.0, 32.0, 26.0, 2.0, 8.0, 36.0, 39.0, 43.0, 41.0, 5.0, 20.0, 17.0, 3.0, 44.0, 43.0, 46.0, 42.0, 14.0, 47.0, 38.0, 31.0, 49.0, 28.0, 15.0, 21.0, 8.0, 14.0]
Len = length(d)

yy = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]

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