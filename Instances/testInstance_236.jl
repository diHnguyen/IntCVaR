edge = [1 2; 1 5; 1 7; 1 9; 1 12; 1 14; 1 15; 2 3; 2 5; 2 8; 2 11; 2 13; 2 14; 3 5; 3 14; 4 5; 4 6; 4 9; 4 11; 4 14; 5 2; 5 3; 6 2; 6 4; 6 10; 6 11; 6 13; 6 15; 7 3; 7 4; 7 5; 7 6; 7 8; 7 12; 7 13; 8 2; 8 3; 8 5; 8 7; 8 12; 8 15; 9 8; 9 13; 10 3; 10 6; 10 11; 10 15; 11 3; 11 5; 11 9; 11 13; 11 14; 11 15; 12 3; 12 5; 12 6; 12 10; 12 11; 12 14; 13 4; 13 9; 13 15; 14 3; 14 9; 14 10; 14 11; 14 13]
cL_orig = [25.0, 181.0, 146.0, 270.0, 368.0, 335.0, 346.0, 45.0, 13.0, 299.0, 65.0, 486.0, 334.0, 75.0, 136.0, 9.0, 52.0, 171.0, 219.0, 438.0, 67.0, 36.0, 134.0, 68.0, 73.0, 59.0, 22.0, 34.0, 101.0, 82.0, 95.0, 3.0, 43.0, 32.0, 289.0, 22.0, 13.0, 30.0, 23.0, 17.0, 4.0, 3.0, 36.0, 114.0, 112.0, 50.0, 11.0, 105.0, 12.0, 87.0, 44.0, 3.0, 53.0, 385.0, 217.0, 138.0, 72.0, 34.0, 0.0, 5.0, 72.0, 29.0, 32.0, 120.0, 149.0, 22.0, 37.0]
cU_orig = [25.0, 181.0, 146.0, 270.0, 368.0, 441.0, 346.0, 45.0, 13.0, 299.0, 65.0, 486.0, 334.0, 75.0, 136.0, 9.0, 92.0, 171.0, 219.0, 438.0, 67.0, 36.0, 134.0, 68.0, 133.0, 59.0, 48.0, 50.0, 101.0, 90.0, 95.0, 3.0, 43.0, 32.0, 289.0, 104.0, 47.0, 30.0, 23.0, 53.0, 4.0, 3.0, 36.0, 492.0, 112.0, 50.0, 27.0, 105.0, 264.0, 87.0, 84.0, 3.0, 53.0, 385.0, 217.0, 138.0, 72.0, 34.0, 80.0, 9.0, 72.0, 29.0, 32.0, 120.0, 149.0, 22.0, 37.0]
d = [9.0, 17.0, 26.0, 31.0, 37.0, 13.0, 6.0, 31.0, 32.0, 37.0, 7.0, 18.0, 21.0, 19.0, 10.0, 19.0, 23.0, 10.0, 24.0, 10.0, 16.0, 37.0, 18.0, 5.0, 36.0, 30.0, 50.0, 40.0, 34.0, 33.0, 11.0, 39.0, 6.0, 49.0, 1.0, 29.0, 21.0, 16.0, 5.0, 39.0, 46.0, 8.0, 46.0, 34.0, 27.0, 15.0, 31.0, 31.0, 23.0, 33.0, 24.0, 15.0, 45.0, 12.0, 12.0, 45.0, 18.0, 27.0, 2.0, 7.0, 22.0, 19.0, 29.0, 38.0, 45.0, 16.0, 31.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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