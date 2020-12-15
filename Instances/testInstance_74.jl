edge = [1 3; 1 5; 1 7; 1 10; 1 11; 1 15; 2 3; 2 6; 2 7; 2 10; 2 11; 3 2; 3 8; 3 10; 3 12; 4 7; 4 11; 4 15; 5 3; 5 4; 5 12; 5 14; 6 3; 6 4; 6 8; 6 10; 6 11; 6 12; 7 2; 7 3; 7 4; 7 10; 7 12; 8 2; 8 3; 8 9; 8 15; 9 4; 9 6; 9 13; 10 3; 10 6; 10 14; 10 15; 11 6; 11 7; 11 10; 11 14; 12 4; 12 7; 12 10; 12 14; 13 4; 13 5; 13 6; 13 7; 13 8; 13 12; 13 15; 14 12]
cL_orig = [49.0, 154.0, 113.0, 225.0, 359.0, 299.0, 38.0, 79.0, 104.0, 336.0, 272.0, 27.0, 188.0, 80.0, 239.0, 27.0, 119.0, 147.0, 85.0, 39.0, 220.0, 426.0, 120.0, 46.0, 24.0, 89.0, 136.0, 297.0, 7.0, 116.0, 10.0, 145.0, 125.0, 167.0, 10.0, 21.0, 331.0, 59.0, 45.0, 44.0, 148.0, 21.0, 16.0, 239.0, 94.0, 156.0, 46.0, 92.0, 268.0, 128.0, 91.0, 25.0, 148.0, 10.0, 327.0, 118.0, 176.0, 44.0, 65.0, 32.0]
cU_orig = [49.0, 154.0, 113.0, 491.0, 359.0, 299.0, 38.0, 79.0, 104.0, 336.0, 332.0, 39.0, 188.0, 584.0, 239.0, 79.0, 119.0, 147.0, 85.0, 39.0, 220.0, 426.0, 120.0, 46.0, 24.0, 177.0, 136.0, 297.0, 7.0, 116.0, 10.0, 145.0, 191.0, 167.0, 474.0, 25.0, 331.0, 59.0, 45.0, 44.0, 520.0, 21.0, 124.0, 239.0, 94.0, 156.0, 46.0, 160.0, 268.0, 128.0, 91.0, 25.0, 148.0, 10.0, 327.0, 118.0, 176.0, 44.0, 113.0, 32.0]
d = [2.0, 47.0, 21.0, 36.0, 24.0, 29.0, 26.0, 23.0, 15.0, 3.0, 32.0, 22.0, 12.0, 18.0, 8.0, 13.0, 15.0, 43.0, 7.0, 9.0, 32.0, 41.0, 22.0, 9.0, 5.0, 1.0, 49.0, 9.0, 18.0, 10.0, 5.0, 5.0, 33.0, 5.0, 28.0, 42.0, 18.0, 49.0, 48.0, 23.0, 27.0, 42.0, 9.0, 14.0, 15.0, 21.0, 45.0, 10.0, 45.0, 21.0, 26.0, 5.0, 20.0, 31.0, 40.0, 16.0, 50.0, 31.0, 31.0, 34.0]
Len = length(d)

yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
