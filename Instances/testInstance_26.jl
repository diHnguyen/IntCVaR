edge = [1 4; 1 5; 1 7; 1 9; 1 15; 2 3; 2 6; 2 7; 2 9; 2 10; 2 14; 3 7; 3 9; 3 11; 3 14; 4 5; 4 7; 4 13; 5 2; 5 8; 5 9; 5 10; 5 14; 5 15; 6 4; 6 5; 6 8; 7 2; 7 9; 7 13; 8 3; 8 4; 8 5; 8 11; 8 13; 9 4; 9 12; 10 3; 10 8; 10 14; 10 15; 11 3; 11 4; 11 8; 11 12; 11 14; 12 2; 12 5; 12 6; 13 7; 13 8; 13 11; 14 2; 14 9; 14 12; 14 13; 14 15]
cL_orig = [110.0, 198.0, 200.0, 281.0, 553.0, 18.0, 140.0, 164.0, 42.0, 296.0, 466.0, 115.0, 114.0, 187.0, 133.0, 4.0, 62.0, 194.0, 11.0, 34.0, 185.0, 56.0, 372.0, 17.0, 40.0, 22.0, 69.0, 47.0, 53.0, 239.0, 45.0, 27.0, 120.0, 19.0, 90.0, 71.0, 33.0, 328.0, 51.0, 120.0, 45.0, 264.0, 154.0, 149.0, 1.0, 7.0, 108.0, 323.0, 169.0, 187.0, 21.0, 57.0, 165.0, 7.0, 65.0, 25.0, 38.0]
cU_orig = [164.0, 198.0, 200.0, 281.0, 553.0, 18.0, 140.0, 164.0, 42.0, 296.0, 466.0, 115.0, 114.0, 187.0, 133.0, 76.0, 62.0, 194.0, 47.0, 34.0, 185.0, 56.0, 426.0, 17.0, 40.0, 22.0, 69.0, 47.0, 53.0, 239.0, 45.0, 205.0, 120.0, 19.0, 90.0, 197.0, 33.0, 328.0, 51.0, 120.0, 45.0, 284.0, 416.0, 149.0, 1.0, 7.0, 108.0, 323.0, 169.0, 187.0, 61.0, 57.0, 165.0, 7.0, 65.0, 25.0, 38.0]
d = [49.0, 39.0, 32.0, 50.0, 32.0, 21.0, 1.0, 34.0, 7.0, 23.0, 50.0, 29.0, 18.0, 49.0, 8.0, 34.0, 19.0, 26.0, 46.0, 47.0, 2.0, 26.0, 27.0, 27.0, 46.0, 32.0, 17.0, 2.0, 14.0, 46.0, 18.0, 29.0, 35.0, 46.0, 41.0, 19.0, 14.0, 37.0, 45.0, 30.0, 12.0, 28.0, 18.0, 22.0, 5.0, 3.0, 36.0, 43.0, 35.0, 11.0, 2.0, 10.0, 14.0, 16.0, 10.0, 12.0, 19.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
