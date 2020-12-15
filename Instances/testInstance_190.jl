edge = [1 4; 1 6; 1 8; 1 10; 1 14; 1 15; 2 3; 2 7; 2 8; 2 11; 2 12; 2 14; 3 9; 3 10; 3 12; 3 14; 3 15; 4 3; 4 6; 5 4; 5 8; 5 9; 5 10; 5 15; 6 2; 6 5; 6 9; 6 11; 6 15; 7 4; 7 5; 7 8; 8 5; 8 6; 8 7; 8 9; 8 10; 8 12; 9 7; 9 14; 10 4; 10 7; 10 8; 10 9; 10 11; 11 5; 11 7; 11 9; 11 10; 11 12; 11 15; 12 8; 12 11; 12 13; 13 8; 13 11; 13 14; 14 5; 14 11; 14 12; 14 15]
cL_orig = [36.0, 235.0, 10.0, 132.0, 212.0, 56.0, 40.0, 24.0, 48.0, 83.0, 209.0, 188.0, 224.0, 136.0, 234.0, 69.0, 579.0, 41.0, 7.0, 48.0, 4.0, 106.0, 87.0, 194.0, 134.0, 47.0, 44.0, 216.0, 39.0, 33.0, 7.0, 45.0, 32.0, 34.0, 35.0, 42.0, 4.0, 33.0, 63.0, 30.0, 219.0, 76.0, 85.0, 27.0, 17.0, 176.0, 89.0, 53.0, 34.0, 12.0, 76.0, 32.0, 33.0, 14.0, 50.0, 19.0, 26.0, 68.0, 62.0, 62.0, 29.0]
cU_orig = [36.0, 235.0, 144.0, 132.0, 574.0, 56.0, 40.0, 24.0, 48.0, 83.0, 209.0, 188.0, 372.0, 136.0, 234.0, 69.0, 579.0, 41.0, 77.0, 48.0, 40.0, 258.0, 87.0, 194.0, 134.0, 47.0, 186.0, 216.0, 39.0, 33.0, 7.0, 45.0, 132.0, 34.0, 35.0, 42.0, 4.0, 209.0, 97.0, 380.0, 219.0, 76.0, 85.0, 27.0, 17.0, 176.0, 89.0, 53.0, 34.0, 12.0, 76.0, 32.0, 45.0, 14.0, 50.0, 19.0, 26.0, 68.0, 62.0, 62.0, 29.0]
d = [10.0, 26.0, 41.0, 40.0, 1.0, 48.0, 42.0, 39.0, 3.0, 18.0, 38.0, 12.0, 8.0, 38.0, 37.0, 48.0, 49.0, 8.0, 7.0, 31.0, 40.0, 47.0, 12.0, 36.0, 22.0, 18.0, 32.0, 34.0, 31.0, 22.0, 29.0, 46.0, 45.0, 47.0, 8.0, 9.0, 46.0, 28.0, 43.0, 25.0, 19.0, 16.0, 2.0, 44.0, 36.0, 48.0, 8.0, 1.0, 13.0, 45.0, 8.0, 9.0, 27.0, 3.0, 17.0, 41.0, 18.0, 19.0, 24.0, 37.0, 50.0]
Len = length(d)

yy = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
