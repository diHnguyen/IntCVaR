edge = [1 5; 1 6; 1 7; 1 10; 1 11; 1 12; 2 4; 2 6; 2 8; 2 11; 2 13; 3 2; 3 4; 3 12; 3 14; 4 7; 4 10; 4 12; 5 7; 5 8; 5 10; 6 2; 6 3; 6 4; 6 9; 6 11; 6 12; 6 13; 6 14; 6 15; 7 5; 7 8; 7 10; 7 12; 8 5; 8 10; 8 11; 9 2; 9 5; 9 6; 9 7; 9 14; 10 2; 10 3; 10 8; 11 9; 11 10; 11 12; 11 14; 11 15; 12 2; 12 4; 12 6; 12 8; 12 9; 13 14; 13 15; 14 7; 14 8; 14 10; 14 12]
cL_orig = [20.0, 198.0, 68.0, 234.0, 128.0, 84.0, 85.0, 51.0, 9.0, 248.0, 404.0, 42.0, 32.0, 83.0, 515.0, 31.0, 125.0, 2.0, 77.0, 28.0, 6.0, 38.0, 80.0, 100.0, 91.0, 89.0, 206.0, 251.0, 104.0, 208.0, 49.0, 5.0, 84.0, 118.0, 143.0, 13.0, 83.0, 334.0, 26.0, 3.0, 22.0, 138.0, 33.0, 75.0, 31.0, 14.0, 1.0, 22.0, 43.0, 2.0, 194.0, 395.0, 10.0, 10.0, 127.0, 16.0, 87.0, 254.0, 147.0, 186.0, 2.0]
cU_orig = [168.0, 198.0, 68.0, 282.0, 128.0, 286.0, 113.0, 51.0, 9.0, 248.0, 404.0, 42.0, 32.0, 83.0, 515.0, 175.0, 125.0, 2.0, 77.0, 28.0, 6.0, 52.0, 80.0, 100.0, 189.0, 91.0, 206.0, 251.0, 104.0, 382.0, 49.0, 5.0, 92.0, 118.0, 143.0, 13.0, 83.0, 334.0, 26.0, 115.0, 22.0, 138.0, 33.0, 75.0, 39.0, 14.0, 7.0, 22.0, 43.0, 28.0, 194.0, 395.0, 44.0, 180.0, 127.0, 16.0, 87.0, 416.0, 147.0, 186.0, 14.0]
d = [50.0, 24.0, 25.0, 20.0, 35.0, 21.0, 17.0, 30.0, 31.0, 5.0, 30.0, 11.0, 42.0, 16.0, 16.0, 35.0, 8.0, 23.0, 12.0, 28.0, 28.0, 7.0, 27.0, 31.0, 27.0, 2.0, 17.0, 30.0, 5.0, 5.0, 5.0, 30.0, 12.0, 18.0, 31.0, 49.0, 12.0, 20.0, 4.0, 16.0, 33.0, 35.0, 46.0, 5.0, 33.0, 38.0, 9.0, 22.0, 4.0, 32.0, 14.0, 20.0, 37.0, 29.0, 21.0, 38.0, 23.0, 31.0, 4.0, 41.0, 32.0]
Len = length(d)

yy = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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