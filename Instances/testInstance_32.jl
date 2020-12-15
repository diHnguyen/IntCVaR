edge = [1 5; 1 6; 1 7; 1 8; 1 13; 1 15; 2 4; 2 5; 2 6; 2 9; 2 10; 2 11; 2 15; 3 4; 3 5; 3 10; 3 15; 4 2; 4 7; 4 14; 5 4; 5 6; 5 10; 6 2; 6 3; 6 10; 6 11; 6 13; 6 14; 6 15; 7 3; 7 5; 7 9; 7 11; 7 13; 7 15; 8 9; 8 10; 8 12; 8 15; 9 3; 9 11; 9 13; 9 14; 10 8; 10 14; 11 6; 11 7; 11 8; 11 15; 12 6; 12 10; 13 6; 13 10; 13 14; 14 2; 14 6; 14 10; 14 13; 14 15]
cL_orig = [87.0, 60.0, 25.0, 106.0, 356.0, 688.0, 36.0, 94.0, 131.0, 98.0, 202.0, 8.0, 512.0, 17.0, 48.0, 223.0, 27.0, 58.0, 12.0, 112.0, 23.0, 50.0, 65.0, 22.0, 2.0, 81.0, 86.0, 115.0, 26.0, 282.0, 55.0, 24.0, 56.0, 54.0, 62.0, 186.0, 27.0, 18.0, 7.0, 164.0, 111.0, 82.0, 20.0, 243.0, 99.0, 179.0, 112.0, 18.0, 137.0, 7.0, 174.0, 43.0, 147.0, 41.0, 0.0, 276.0, 173.0, 62.0, 38.0, 16.0]
cU_orig = [181.0, 60.0, 25.0, 106.0, 356.0, 688.0, 36.0, 94.0, 131.0, 130.0, 202.0, 8.0, 512.0, 37.0, 48.0, 245.0, 65.0, 58.0, 12.0, 112.0, 23.0, 50.0, 65.0, 22.0, 14.0, 235.0, 86.0, 157.0, 62.0, 568.0, 113.0, 24.0, 56.0, 74.0, 62.0, 456.0, 27.0, 18.0, 7.0, 164.0, 267.0, 82.0, 20.0, 243.0, 99.0, 179.0, 112.0, 62.0, 157.0, 7.0, 174.0, 43.0, 147.0, 41.0, 4.0, 276.0, 173.0, 166.0, 58.0, 28.0]
d = [45.0, 31.0, 46.0, 42.0, 28.0, 13.0, 8.0, 5.0, 22.0, 13.0, 44.0, 43.0, 46.0, 50.0, 4.0, 22.0, 5.0, 2.0, 46.0, 50.0, 2.0, 27.0, 14.0, 34.0, 31.0, 11.0, 32.0, 43.0, 27.0, 32.0, 16.0, 10.0, 31.0, 22.0, 46.0, 35.0, 28.0, 32.0, 30.0, 12.0, 4.0, 7.0, 6.0, 44.0, 11.0, 15.0, 7.0, 20.0, 39.0, 13.0, 36.0, 29.0, 40.0, 50.0, 41.0, 15.0, 3.0, 19.0, 2.0, 47.0]
Len = length(d)

yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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