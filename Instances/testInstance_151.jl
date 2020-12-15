edge = [1 6; 1 7; 1 9; 1 10; 1 15; 2 3; 2 4; 2 6; 2 10; 2 11; 2 12; 3 5; 3 7; 4 2; 4 9; 4 12; 4 14; 5 4; 5 11; 5 15; 6 3; 6 4; 6 5; 6 8; 6 9; 6 11; 6 15; 7 10; 7 11; 8 5; 8 9; 8 13; 8 15; 9 3; 9 6; 9 7; 9 8; 9 10; 10 2; 10 3; 10 11; 10 12; 10 15; 11 3; 11 4; 11 5; 11 6; 11 9; 11 10; 11 12; 12 2; 12 3; 12 6; 12 7; 12 13; 12 14; 12 15; 13 15; 14 4; 14 6; 14 15]
cL_orig = [172.0, 134.0, 17.0, 89.0, 257.0, 46.0, 33.0, 135.0, 10.0, 168.0, 396.0, 61.0, 87.0, 6.0, 211.0, 115.0, 149.0, 27.0, 44.0, 348.0, 12.0, 26.0, 31.0, 98.0, 60.0, 56.0, 304.0, 119.0, 4.0, 114.0, 12.0, 133.0, 151.0, 284.0, 9.0, 27.0, 14.0, 34.0, 169.0, 309.0, 8.0, 69.0, 65.0, 348.0, 22.0, 18.0, 245.0, 27.0, 10.0, 48.0, 353.0, 187.0, 122.0, 90.0, 35.0, 2.0, 39.0, 66.0, 67.0, 157.0, 27.0]
cU_orig = [172.0, 134.0, 45.0, 111.0, 257.0, 46.0, 143.0, 135.0, 10.0, 168.0, 396.0, 61.0, 87.0, 6.0, 211.0, 115.0, 605.0, 27.0, 44.0, 628.0, 12.0, 34.0, 31.0, 98.0, 76.0, 264.0, 304.0, 119.0, 176.0, 114.0, 12.0, 133.0, 151.0, 284.0, 9.0, 27.0, 14.0, 44.0, 169.0, 309.0, 8.0, 69.0, 79.0, 348.0, 486.0, 18.0, 245.0, 151.0, 10.0, 48.0, 353.0, 187.0, 314.0, 90.0, 55.0, 6.0, 39.0, 66.0, 265.0, 157.0, 27.0]
d = [49.0, 17.0, 30.0, 26.0, 46.0, 29.0, 12.0, 28.0, 6.0, 21.0, 50.0, 33.0, 9.0, 37.0, 22.0, 42.0, 2.0, 19.0, 10.0, 24.0, 39.0, 43.0, 2.0, 6.0, 35.0, 11.0, 33.0, 13.0, 15.0, 7.0, 49.0, 42.0, 44.0, 7.0, 49.0, 43.0, 34.0, 42.0, 41.0, 10.0, 40.0, 12.0, 20.0, 31.0, 44.0, 23.0, 28.0, 41.0, 48.0, 47.0, 4.0, 7.0, 41.0, 14.0, 42.0, 23.0, 48.0, 31.0, 38.0, 23.0, 47.0]
Len = length(d)

yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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