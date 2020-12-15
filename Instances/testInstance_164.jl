edge = [1 2; 1 4; 1 5; 1 7; 1 10; 2 3; 2 7; 2 9; 2 12; 3 4; 3 5; 3 7; 3 10; 3 12; 4 7; 4 11; 5 2; 5 10; 5 12; 5 13; 6 2; 6 3; 6 13; 7 2; 7 3; 7 6; 7 10; 8 2; 8 4; 8 10; 8 11; 8 12; 9 2; 9 4; 9 7; 9 13; 9 15; 10 4; 10 8; 10 13; 11 5; 11 10; 11 13; 12 2; 12 8; 12 9; 13 4; 13 5; 13 9; 14 3; 14 5; 14 7; 14 8]
cL_orig = [37.0, 110.0, 48.0, 29.0, 182.0, 39.0, 33.0, 111.0, 387.0, 33.0, 23.0, 79.0, 0.0, 386.0, 115.0, 70.0, 82.0, 24.0, 245.0, 48.0, 145.0, 120.0, 265.0, 115.0, 182.0, 25.0, 4.0, 269.0, 2.0, 2.0, 71.0, 152.0, 62.0, 211.0, 68.0, 90.0, 263.0, 126.0, 29.0, 120.0, 269.0, 4.0, 24.0, 66.0, 132.0, 55.0, 33.0, 130.0, 115.0, 39.0, 4.0, 350.0, 272.0]
cU_orig = [37.0, 110.0, 48.0, 29.0, 182.0, 39.0, 33.0, 111.0, 387.0, 33.0, 23.0, 79.0, 2.0, 386.0, 115.0, 70.0, 82.0, 24.0, 245.0, 48.0, 145.0, 120.0, 265.0, 115.0, 182.0, 25.0, 4.0, 269.0, 222.0, 158.0, 219.0, 152.0, 374.0, 211.0, 68.0, 90.0, 263.0, 126.0, 29.0, 120.0, 269.0, 88.0, 24.0, 74.0, 132.0, 55.0, 33.0, 178.0, 115.0, 93.0, 22.0, 350.0, 272.0]
d = [40.0, 41.0, 41.0, 38.0, 4.0, 11.0, 34.0, 5.0, 48.0, 46.0, 24.0, 30.0, 10.0, 14.0, 38.0, 19.0, 22.0, 36.0, 29.0, 24.0, 24.0, 13.0, 43.0, 30.0, 32.0, 16.0, 29.0, 32.0, 31.0, 3.0, 41.0, 14.0, 14.0, 1.0, 47.0, 18.0, 10.0, 30.0, 47.0, 24.0, 50.0, 47.0, 45.0, 48.0, 14.0, 15.0, 50.0, 23.0, 1.0, 26.0, 38.0, 46.0, 19.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
