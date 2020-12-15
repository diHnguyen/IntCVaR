edge = [1 3; 1 9; 1 11; 1 12; 1 14; 2 3; 2 5; 2 8; 2 9; 2 10; 2 11; 2 14; 2 15; 3 10; 3 11; 3 13; 4 2; 4 7; 4 8; 5 11; 5 12; 5 15; 6 5; 6 7; 6 8; 6 9; 6 10; 6 11; 6 15; 7 2; 7 6; 7 8; 7 13; 7 15; 8 3; 8 5; 8 12; 8 13; 8 14; 9 8; 9 10; 9 13; 9 14; 10 7; 10 8; 10 9; 11 3; 11 8; 11 12; 12 5; 12 8; 12 9; 12 10; 13 6; 13 9; 13 14; 13 15; 14 6; 14 8; 14 10; 14 11; 14 12; 14 13; 14 15]
cL_orig = [81.0, 130.0, 288.0, 177.0, 574.0, 28.0, 119.0, 250.0, 135.0, 18.0, 212.0, 292.0, 351.0, 13.0, 215.0, 22.0, 90.0, 58.0, 136.0, 13.0, 79.0, 227.0, 0.0, 16.0, 79.0, 136.0, 90.0, 10.0, 298.0, 99.0, 0.0, 43.0, 25.0, 243.0, 122.0, 112.0, 86.0, 114.0, 67.0, 12.0, 13.0, 19.0, 11.0, 19.0, 23.0, 25.0, 315.0, 129.0, 10.0, 41.0, 184.0, 32.0, 67.0, 11.0, 30.0, 43.0, 68.0, 254.0, 46.0, 54.0, 34.0, 21.0, 39.0, 10.0]
cU_orig = [81.0, 130.0, 570.0, 177.0, 574.0, 38.0, 119.0, 250.0, 135.0, 254.0, 212.0, 292.0, 787.0, 13.0, 215.0, 304.0, 90.0, 164.0, 172.0, 13.0, 79.0, 227.0, 52.0, 16.0, 79.0, 136.0, 108.0, 86.0, 298.0, 195.0, 54.0, 43.0, 261.0, 243.0, 182.0, 112.0, 86.0, 114.0, 67.0, 12.0, 13.0, 49.0, 43.0, 19.0, 23.0, 25.0, 315.0, 129.0, 10.0, 41.0, 184.0, 32.0, 67.0, 11.0, 30.0, 43.0, 68.0, 254.0, 80.0, 82.0, 34.0, 67.0, 39.0, 10.0]
d = [17.0, 14.0, 38.0, 42.0, 49.0, 47.0, 19.0, 14.0, 16.0, 46.0, 7.0, 34.0, 34.0, 23.0, 43.0, 50.0, 36.0, 21.0, 13.0, 22.0, 4.0, 35.0, 50.0, 25.0, 49.0, 34.0, 23.0, 36.0, 25.0, 17.0, 43.0, 40.0, 6.0, 11.0, 10.0, 28.0, 7.0, 31.0, 18.0, 38.0, 48.0, 43.0, 45.0, 17.0, 14.0, 28.0, 12.0, 42.0, 22.0, 38.0, 29.0, 46.0, 22.0, 20.0, 46.0, 29.0, 13.0, 23.0, 17.0, 22.0, 4.0, 18.0, 6.0, 46.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]

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
