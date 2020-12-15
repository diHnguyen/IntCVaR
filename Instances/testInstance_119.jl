edge = [1 3; 1 9; 1 10; 1 14; 2 3; 2 8; 2 11; 2 12; 2 15; 3 4; 3 8; 3 9; 3 10; 3 12; 3 14; 4 2; 4 7; 4 10; 5 2; 5 3; 5 9; 5 11; 6 2; 7 3; 7 5; 7 6; 7 12; 7 14; 8 2; 8 4; 8 7; 9 2; 9 4; 9 12; 9 13; 10 4; 10 8; 10 14; 10 15; 11 8; 11 9; 11 13; 12 6; 12 8; 12 11; 12 13; 13 5; 13 7; 13 8; 13 11; 13 14; 13 15; 14 6]
cL_orig = [90.0, 84.0, 41.0, 9.0, 43.0, 172.0, 164.0, 64.0, 488.0, 28.0, 21.0, 47.0, 0.0, 111.0, 375.0, 12.0, 58.0, 6.0, 111.0, 41.0, 98.0, 278.0, 134.0, 82.0, 46.0, 21.0, 147.0, 35.0, 202.0, 146.0, 15.0, 164.0, 6.0, 20.0, 77.0, 171.0, 2.0, 127.0, 54.0, 137.0, 93.0, 61.0, 120.0, 42.0, 39.0, 49.0, 274.0, 3.0, 159.0, 44.0, 1.0, 18.0, 293.0]
cU_orig = [90.0, 84.0, 41.0, 9.0, 43.0, 172.0, 164.0, 64.0, 488.0, 28.0, 23.0, 257.0, 24.0, 385.0, 375.0, 24.0, 58.0, 6.0, 111.0, 41.0, 294.0, 312.0, 134.0, 82.0, 46.0, 21.0, 353.0, 35.0, 202.0, 146.0, 15.0, 164.0, 246.0, 120.0, 77.0, 171.0, 2.0, 127.0, 54.0, 137.0, 93.0, 61.0, 120.0, 42.0, 39.0, 49.0, 312.0, 3.0, 159.0, 44.0, 1.0, 18.0, 293.0]
d = [5.0, 32.0, 42.0, 1.0, 38.0, 18.0, 34.0, 7.0, 39.0, 19.0, 45.0, 44.0, 12.0, 23.0, 26.0, 10.0, 44.0, 32.0, 27.0, 46.0, 13.0, 50.0, 38.0, 37.0, 28.0, 3.0, 34.0, 29.0, 29.0, 16.0, 6.0, 5.0, 9.0, 37.0, 11.0, 12.0, 1.0, 28.0, 28.0, 32.0, 47.0, 18.0, 47.0, 30.0, 16.0, 48.0, 19.0, 35.0, 41.0, 20.0, 28.0, 50.0, 4.0]
Len = length(d)

yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
