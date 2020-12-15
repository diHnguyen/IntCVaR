edge = [1 2; 1 3; 1 5; 1 6; 1 7; 1 10; 1 13; 2 3; 2 4; 2 5; 2 8; 2 10; 2 11; 2 13; 3 2; 3 11; 3 14; 4 8; 4 9; 4 11; 4 12; 5 3; 5 8; 5 11; 5 12; 6 2; 6 5; 6 7; 6 8; 6 9; 6 14; 7 3; 7 9; 7 11; 7 13; 7 14; 7 15; 8 2; 8 6; 8 10; 8 12; 9 2; 9 3; 9 4; 9 7; 10 5; 10 15; 11 2; 11 4; 11 6; 11 7; 11 12; 11 14; 11 15; 12 5; 12 7; 12 8; 12 14; 13 2; 13 4; 13 5; 13 8; 13 10; 13 12; 14 2; 14 7; 14 10; 14 11]
cL_orig = [31.0, 83.0, 35.0, 233.0, 241.0, 394.0, 301.0, 4.0, 76.0, 54.0, 78.0, 33.0, 69.0, 515.0, 26.0, 292.0, 425.0, 88.0, 55.0, 177.0, 377.0, 53.0, 45.0, 248.0, 3.0, 29.0, 0.0, 10.0, 85.0, 128.0, 338.0, 141.0, 34.0, 23.0, 126.0, 259.0, 67.0, 240.0, 90.0, 63.0, 51.0, 125.0, 104.0, 159.0, 65.0, 239.0, 228.0, 243.0, 317.0, 43.0, 43.0, 22.0, 20.0, 27.0, 37.0, 202.0, 101.0, 1.0, 193.0, 124.0, 135.0, 178.0, 110.0, 2.0, 38.0, 35.0, 168.0, 21.0]
cU_orig = [31.0, 83.0, 35.0, 233.0, 241.0, 394.0, 301.0, 4.0, 76.0, 54.0, 78.0, 33.0, 711.0, 515.0, 26.0, 292.0, 425.0, 88.0, 215.0, 177.0, 377.0, 113.0, 45.0, 248.0, 3.0, 29.0, 50.0, 78.0, 85.0, 128.0, 338.0, 141.0, 34.0, 23.0, 126.0, 259.0, 553.0, 240.0, 90.0, 63.0, 77.0, 125.0, 104.0, 159.0, 65.0, 239.0, 228.0, 243.0, 317.0, 43.0, 43.0, 22.0, 74.0, 189.0, 89.0, 232.0, 101.0, 1.0, 193.0, 124.0, 183.0, 318.0, 110.0, 2.0, 90.0, 35.0, 168.0, 25.0]
d = [50.0, 30.0, 40.0, 12.0, 41.0, 24.0, 44.0, 29.0, 21.0, 43.0, 46.0, 41.0, 45.0, 25.0, 2.0, 2.0, 43.0, 24.0, 5.0, 8.0, 9.0, 48.0, 49.0, 16.0, 22.0, 34.0, 30.0, 48.0, 31.0, 45.0, 13.0, 1.0, 40.0, 37.0, 1.0, 2.0, 12.0, 40.0, 32.0, 4.0, 11.0, 29.0, 33.0, 12.0, 30.0, 10.0, 38.0, 29.0, 30.0, 27.0, 6.0, 47.0, 9.0, 16.0, 11.0, 12.0, 2.0, 48.0, 26.0, 6.0, 26.0, 44.0, 29.0, 43.0, 22.0, 22.0, 35.0, 9.0]
Len = length(d)

yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]

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