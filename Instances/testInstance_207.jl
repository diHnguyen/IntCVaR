edge = [1 4; 1 5; 1 7; 1 9; 1 13; 2 7; 2 11; 2 13; 2 15; 3 6; 3 13; 3 15; 4 2; 4 3; 4 10; 4 14; 5 2; 5 4; 5 6; 5 7; 5 10; 5 14; 5 15; 6 3; 6 8; 6 9; 6 14; 7 5; 7 6; 7 9; 7 11; 8 4; 8 7; 8 11; 8 12; 9 2; 9 5; 9 7; 9 14; 9 15; 10 4; 10 5; 10 8; 10 11; 10 13; 10 15; 11 4; 11 5; 11 6; 11 9; 11 10; 11 12; 11 13; 11 14; 12 2; 12 4; 12 6; 12 8; 12 9; 12 14; 13 2; 13 3; 13 9; 13 10; 14 2; 14 4; 14 11; 14 13]
cL_orig = [20.0, 22.0, 59.0, 0.0, 124.0, 176.0, 23.0, 338.0, 29.0, 31.0, 283.0, 416.0, 59.0, 34.0, 54.0, 329.0, 102.0, 38.0, 9.0, 88.0, 2.0, 236.0, 75.0, 57.0, 4.0, 87.0, 246.0, 14.0, 18.0, 99.0, 19.0, 45.0, 7.0, 123.0, 15.0, 238.0, 13.0, 79.0, 219.0, 257.0, 124.0, 124.0, 0.0, 43.0, 32.0, 91.0, 239.0, 299.0, 11.0, 19.0, 11.0, 33.0, 51.0, 69.0, 284.0, 220.0, 16.0, 96.0, 72.0, 57.0, 64.0, 411.0, 83.0, 5.0, 17.0, 89.0, 46.0, 19.0]
cU_orig = [20.0, 266.0, 59.0, 488.0, 124.0, 176.0, 155.0, 338.0, 83.0, 31.0, 283.0, 416.0, 59.0, 34.0, 100.0, 329.0, 102.0, 38.0, 19.0, 88.0, 2.0, 352.0, 219.0, 57.0, 16.0, 117.0, 246.0, 14.0, 18.0, 99.0, 19.0, 131.0, 7.0, 123.0, 15.0, 238.0, 13.0, 79.0, 219.0, 257.0, 124.0, 124.0, 136.0, 43.0, 210.0, 91.0, 239.0, 299.0, 11.0, 19.0, 11.0, 33.0, 59.0, 69.0, 284.0, 220.0, 16.0, 200.0, 72.0, 57.0, 472.0, 411.0, 83.0, 5.0, 129.0, 89.0, 254.0, 19.0]
d = [18.0, 26.0, 15.0, 6.0, 16.0, 41.0, 3.0, 11.0, 6.0, 34.0, 3.0, 29.0, 8.0, 50.0, 39.0, 43.0, 12.0, 42.0, 38.0, 38.0, 34.0, 13.0, 17.0, 39.0, 46.0, 21.0, 10.0, 8.0, 31.0, 22.0, 29.0, 12.0, 44.0, 48.0, 48.0, 16.0, 16.0, 9.0, 22.0, 28.0, 50.0, 15.0, 50.0, 49.0, 13.0, 49.0, 32.0, 3.0, 29.0, 7.0, 46.0, 16.0, 14.0, 7.0, 47.0, 34.0, 7.0, 2.0, 16.0, 44.0, 33.0, 8.0, 1.0, 3.0, 34.0, 25.0, 38.0, 46.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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