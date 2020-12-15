edge = [1 2; 1 4; 1 5; 1 7; 1 11; 1 15; 2 4; 2 6; 2 7; 2 8; 2 9; 2 10; 2 11; 2 15; 3 2; 3 8; 3 13; 4 7; 4 8; 4 10; 4 11; 4 12; 4 15; 5 2; 5 3; 5 6; 5 11; 5 12; 6 8; 6 11; 6 15; 7 2; 7 4; 7 5; 7 9; 7 13; 8 3; 8 4; 8 7; 8 13; 8 15; 9 2; 9 4; 9 5; 9 7; 9 11; 10 3; 10 4; 10 6; 10 13; 10 14; 10 15; 11 2; 11 9; 12 7; 12 9; 12 11; 12 14; 13 11; 14 2; 14 3; 14 8; 14 10; 14 12; 14 13]
cL_orig = [15.0, 123.0, 47.0, 42.0, 76.0, 222.0, 18.0, 121.0, 119.0, 188.0, 193.0, 6.0, 330.0, 639.0, 11.0, 31.0, 289.0, 0.0, 113.0, 0.0, 231.0, 388.0, 325.0, 65.0, 2.0, 12.0, 74.0, 73.0, 1.0, 1.0, 167.0, 192.0, 56.0, 68.0, 2.0, 72.0, 123.0, 139.0, 13.0, 19.0, 35.0, 62.0, 42.0, 107.0, 63.0, 65.0, 173.0, 266.0, 144.0, 42.0, 62.0, 71.0, 220.0, 53.0, 150.0, 22.0, 8.0, 0.0, 2.0, 107.0, 311.0, 61.0, 191.0, 26.0, 6.0]
cU_orig = [15.0, 123.0, 47.0, 58.0, 314.0, 222.0, 26.0, 121.0, 119.0, 188.0, 193.0, 6.0, 330.0, 639.0, 29.0, 121.0, 289.0, 6.0, 113.0, 36.0, 419.0, 388.0, 325.0, 65.0, 4.0, 12.0, 74.0, 73.0, 105.0, 1.0, 167.0, 192.0, 170.0, 68.0, 10.0, 342.0, 123.0, 139.0, 13.0, 369.0, 35.0, 62.0, 66.0, 107.0, 63.0, 99.0, 189.0, 266.0, 144.0, 42.0, 62.0, 71.0, 220.0, 91.0, 150.0, 22.0, 8.0, 76.0, 6.0, 107.0, 311.0, 61.0, 191.0, 26.0, 6.0]
d = [14.0, 40.0, 8.0, 21.0, 10.0, 17.0, 2.0, 43.0, 33.0, 5.0, 1.0, 16.0, 7.0, 44.0, 19.0, 20.0, 30.0, 37.0, 28.0, 21.0, 18.0, 17.0, 24.0, 23.0, 45.0, 49.0, 3.0, 19.0, 3.0, 38.0, 32.0, 36.0, 8.0, 34.0, 2.0, 18.0, 19.0, 4.0, 4.0, 38.0, 49.0, 18.0, 14.0, 21.0, 47.0, 46.0, 31.0, 29.0, 48.0, 26.0, 18.0, 46.0, 23.0, 50.0, 10.0, 49.0, 33.0, 27.0, 8.0, 29.0, 24.0, 37.0, 49.0, 34.0, 22.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
