edge = [1 2; 1 3; 1 6; 1 7; 1 8; 1 13; 1 15; 2 6; 2 9; 2 11; 2 12; 2 14; 3 4; 3 8; 3 9; 3 12; 3 15; 4 3; 4 6; 4 7; 4 8; 4 12; 4 13; 5 12; 5 14; 5 15; 6 5; 6 7; 6 15; 7 6; 7 9; 7 10; 7 12; 7 15; 8 2; 8 5; 8 11; 8 13; 9 3; 9 10; 9 12; 9 13; 9 15; 10 2; 10 3; 10 4; 10 5; 10 6; 10 11; 10 14; 11 3; 11 12; 11 15; 12 2; 12 3; 12 5; 12 8; 12 11; 12 14; 12 15; 13 2; 13 10; 13 11; 14 2; 14 5; 14 6; 14 8; 14 10; 14 11; 14 12]
cL_orig = [43.0, 49.0, 91.0, 149.0, 49.0, 89.0, 283.0, 95.0, 59.0, 346.0, 265.0, 169.0, 24.0, 56.0, 6.0, 436.0, 163.0, 11.0, 44.0, 147.0, 10.0, 72.0, 411.0, 60.0, 9.0, 252.0, 30.0, 1.0, 292.0, 37.0, 22.0, 127.0, 53.0, 1.0, 137.0, 58.0, 26.0, 5.0, 109.0, 12.0, 30.0, 139.0, 128.0, 136.0, 10.0, 46.0, 33.0, 42.0, 7.0, 145.0, 90.0, 13.0, 53.0, 326.0, 285.0, 123.0, 197.0, 3.0, 4.0, 33.0, 523.0, 99.0, 49.0, 35.0, 339.0, 43.0, 55.0, 104.0, 77.0, 84.0]
cU_orig = [43.0, 49.0, 91.0, 149.0, 541.0, 89.0, 771.0, 95.0, 587.0, 346.0, 265.0, 169.0, 24.0, 56.0, 6.0, 436.0, 163.0, 89.0, 44.0, 147.0, 10.0, 72.0, 411.0, 326.0, 25.0, 252.0, 30.0, 37.0, 292.0, 37.0, 22.0, 127.0, 437.0, 211.0, 137.0, 58.0, 64.0, 201.0, 109.0, 12.0, 114.0, 173.0, 132.0, 136.0, 42.0, 46.0, 33.0, 42.0, 7.0, 145.0, 90.0, 13.0, 53.0, 326.0, 285.0, 479.0, 197.0, 3.0, 184.0, 213.0, 523.0, 99.0, 49.0, 35.0, 339.0, 43.0, 193.0, 104.0, 77.0, 84.0]
d = [15.0, 48.0, 34.0, 50.0, 8.0, 5.0, 8.0, 15.0, 17.0, 7.0, 44.0, 24.0, 17.0, 9.0, 10.0, 2.0, 6.0, 22.0, 15.0, 38.0, 25.0, 18.0, 33.0, 4.0, 8.0, 18.0, 46.0, 26.0, 42.0, 30.0, 50.0, 31.0, 17.0, 19.0, 13.0, 34.0, 42.0, 27.0, 22.0, 2.0, 25.0, 50.0, 39.0, 10.0, 21.0, 19.0, 14.0, 22.0, 3.0, 20.0, 11.0, 39.0, 23.0, 24.0, 20.0, 40.0, 13.0, 18.0, 37.0, 48.0, 16.0, 38.0, 3.0, 46.0, 16.0, 40.0, 23.0, 42.0, 36.0, 48.0]
Len = length(d)

yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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