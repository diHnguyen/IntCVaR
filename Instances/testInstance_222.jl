edge = [1 3; 1 10; 1 12; 1 15; 2 4; 2 6; 2 7; 2 9; 2 12; 2 15; 3 6; 3 8; 3 9; 3 11; 4 3; 4 6; 4 7; 4 9; 4 13; 4 14; 5 7; 5 9; 6 4; 6 5; 6 8; 6 9; 7 3; 7 4; 7 8; 7 10; 7 13; 7 15; 9 10; 9 13; 10 2; 10 3; 10 4; 10 5; 10 8; 10 15; 11 5; 11 6; 11 8; 11 13; 12 3; 12 8; 12 14; 13 3; 13 4; 13 6; 13 9; 13 12; 13 15; 14 3; 14 4; 14 8; 14 9; 14 10]
cL_orig = [57.0, 5.0, 139.0, 50.0, 17.0, 6.0, 217.0, 0.0, 375.0, 23.0, 91.0, 2.0, 127.0, 239.0, 40.0, 44.0, 116.0, 128.0, 13.0, 122.0, 75.0, 1.0, 90.0, 22.0, 64.0, 26.0, 8.0, 41.0, 3.0, 60.0, 0.0, 80.0, 9.0, 83.0, 214.0, 104.0, 0.0, 85.0, 55.0, 147.0, 69.0, 91.0, 43.0, 23.0, 372.0, 121.0, 51.0, 229.0, 147.0, 15.0, 118.0, 44.0, 57.0, 111.0, 46.0, 233.0, 33.0, 141.0]
cU_orig = [57.0, 5.0, 721.0, 50.0, 17.0, 84.0, 217.0, 2.0, 375.0, 1195.0, 91.0, 452.0, 127.0, 339.0, 40.0, 48.0, 116.0, 128.0, 45.0, 522.0, 75.0, 7.0, 90.0, 22.0, 64.0, 26.0, 14.0, 41.0, 3.0, 60.0, 102.0, 516.0, 23.0, 83.0, 214.0, 104.0, 306.0, 175.0, 55.0, 147.0, 69.0, 139.0, 43.0, 23.0, 372.0, 121.0, 97.0, 229.0, 147.0, 41.0, 118.0, 44.0, 57.0, 111.0, 46.0, 233.0, 45.0, 199.0]
d = [44.0, 27.0, 35.0, 13.0, 21.0, 46.0, 15.0, 38.0, 34.0, 41.0, 14.0, 48.0, 34.0, 37.0, 34.0, 21.0, 17.0, 47.0, 19.0, 37.0, 35.0, 22.0, 42.0, 7.0, 39.0, 23.0, 17.0, 8.0, 29.0, 7.0, 24.0, 17.0, 17.0, 9.0, 38.0, 37.0, 29.0, 23.0, 46.0, 20.0, 28.0, 40.0, 19.0, 14.0, 45.0, 22.0, 25.0, 35.0, 37.0, 21.0, 48.0, 44.0, 28.0, 42.0, 23.0, 11.0, 11.0, 14.0]
Len = length(d)

yy = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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