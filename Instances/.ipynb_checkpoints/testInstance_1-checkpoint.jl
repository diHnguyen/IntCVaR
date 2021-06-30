edge = [1 3; 1 4; 1 5; 1 6; 1 10; 2 5; 2 6; 2 7; 2 10; 2 12; 2 13; 2 15; 3 9; 3 10; 3 15; 4 3; 4 14; 5 7; 5 8; 5 9; 5 10; 5 12; 5 13; 5 15; 6 2; 6 7; 6 8; 6 9; 6 10; 7 3; 7 6; 7 8; 8 2; 8 3; 8 4; 8 5; 8 12; 8 13; 9 2; 9 4; 9 5; 9 8; 9 11; 9 12; 9 13; 10 4; 10 8; 10 13; 11 3; 11 6; 11 9; 11 12; 12 3; 12 5; 12 10; 13 2; 13 4; 13 5; 13 6; 13 12; 13 14; 13 15; 14 4; 14 5; 14 7; 14 8; 14 11; 14 12; 14 15]
cL_orig = [56.0, 51.0, 182.0, 73.0, 13.0, 29.0, 42.0, 250.0, 9.0, 316.0, 95.0, 177.0, 21.0, 52.0, 318.0, 33.0, 415.0, 26.0, 37.0, 1.0, 17.0, 128.0, 356.0, 407.0, 23.0, 44.0, 100.0, 49.0, 172.0, 6.0, 21.0, 8.0, 44.0, 47.0, 32.0, 41.0, 99.0, 29.0, 284.0, 157.0, 5.0, 44.0, 32.0, 13.0, 47.0, 143.0, 61.0, 58.0, 391.0, 136.0, 54.0, 24.0, 164.0, 118.0, 93.0, 211.0, 295.0, 65.0, 228.0, 48.0, 14.0, 69.0, 193.0, 44.0, 72.0, 36.0, 79.0, 14.0, 42.0]
cU_orig = [56.0, 63.0, 182.0, 187.0, 45.0, 29.0, 56.0, 250.0, 327.0, 316.0, 95.0, 177.0, 21.0, 52.0, 318.0, 33.0, 415.0, 34.0, 111.0, 211.0, 17.0, 202.0, 356.0, 407.0, 167.0, 44.0, 100.0, 49.0, 172.0, 32.0, 21.0, 8.0, 44.0, 47.0, 32.0, 225.0, 99.0, 29.0, 284.0, 157.0, 143.0, 44.0, 84.0, 13.0, 93.0, 143.0, 61.0, 58.0, 391.0, 136.0, 54.0, 68.0, 164.0, 196.0, 93.0, 215.0, 295.0, 525.0, 228.0, 48.0, 14.0, 93.0, 193.0, 44.0, 72.0, 36.0, 107.0, 14.0, 42.0]
d = [3.0, 30.0, 31.0, 13.0, 23.0, 17.0, 9.0, 26.0, 32.0, 33.0, 13.0, 29.0, 13.0, 33.0, 39.0, 32.0, 41.0, 27.0, 41.0, 9.0, 44.0, 2.0, 46.0, 13.0, 23.0, 22.0, 33.0, 6.0, 37.0, 21.0, 13.0, 24.0, 48.0, 44.0, 14.0, 37.0, 30.0, 38.0, 46.0, 13.0, 50.0, 37.0, 25.0, 1.0, 45.0, 27.0, 40.0, 31.0, 40.0, 11.0, 43.0, 17.0, 49.0, 17.0, 19.0, 21.0, 2.0, 33.0, 12.0, 30.0, 19.0, 2.0, 33.0, 29.0, 36.0, 8.0, 46.0, 45.0, 22.0]
Len = length(d)

yy = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1]

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