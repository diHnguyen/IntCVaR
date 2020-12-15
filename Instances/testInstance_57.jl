edge = [1 6; 1 7; 1 8; 1 13; 2 4; 2 6; 2 8; 2 11; 2 12; 3 2; 3 5; 3 7; 3 9; 3 12; 3 13; 4 8; 4 10; 4 12; 4 14; 5 4; 5 6; 5 7; 5 8; 5 9; 5 11; 5 12; 6 2; 6 7; 6 9; 6 14; 7 3; 7 5; 7 6; 7 10; 7 11; 7 12; 7 13; 8 5; 8 9; 8 10; 8 11; 8 13; 8 14; 8 15; 9 3; 9 4; 9 6; 9 7; 9 8; 9 10; 9 14; 10 2; 10 4; 10 11; 11 2; 11 6; 11 7; 11 12; 12 7; 12 9; 12 10; 12 11; 13 3; 13 9; 13 10; 13 12; 14 6; 14 8; 14 11; 14 15]
cL_orig = [127.0, 246.0, 227.0, 364.0, 81.0, 20.0, 108.0, 106.0, 281.0, 38.0, 39.0, 38.0, 171.0, 175.0, 220.0, 56.0, 196.0, 242.0, 35.0, 21.0, 3.0, 49.0, 82.0, 196.0, 151.0, 131.0, 96.0, 20.0, 44.0, 130.0, 21.0, 52.0, 19.0, 94.0, 94.0, 186.0, 62.0, 61.0, 44.0, 0.0, 73.0, 23.0, 128.0, 232.0, 75.0, 3.0, 13.0, 55.0, 29.0, 50.0, 28.0, 200.0, 23.0, 15.0, 101.0, 89.0, 80.0, 5.0, 130.0, 106.0, 89.0, 19.0, 32.0, 166.0, 34.0, 6.0, 20.0, 24.0, 124.0, 16.0]
cU_orig = [127.0, 246.0, 227.0, 364.0, 81.0, 20.0, 108.0, 726.0, 281.0, 38.0, 39.0, 38.0, 171.0, 175.0, 220.0, 56.0, 196.0, 242.0, 35.0, 21.0, 5.0, 49.0, 82.0, 196.0, 151.0, 199.0, 96.0, 20.0, 44.0, 130.0, 21.0, 52.0, 57.0, 118.0, 94.0, 284.0, 62.0, 223.0, 44.0, 20.0, 73.0, 23.0, 128.0, 232.0, 75.0, 53.0, 19.0, 121.0, 35.0, 50.0, 58.0, 200.0, 57.0, 33.0, 665.0, 89.0, 80.0, 61.0, 368.0, 106.0, 89.0, 19.0, 32.0, 166.0, 34.0, 58.0, 20.0, 24.0, 124.0, 18.0]
d = [33.0, 25.0, 6.0, 39.0, 39.0, 10.0, 8.0, 39.0, 25.0, 42.0, 4.0, 24.0, 30.0, 35.0, 2.0, 45.0, 5.0, 49.0, 46.0, 22.0, 37.0, 23.0, 29.0, 40.0, 9.0, 16.0, 1.0, 47.0, 36.0, 37.0, 10.0, 27.0, 32.0, 43.0, 8.0, 46.0, 40.0, 28.0, 11.0, 25.0, 18.0, 12.0, 21.0, 23.0, 42.0, 9.0, 50.0, 17.0, 21.0, 14.0, 24.0, 2.0, 47.0, 1.0, 37.0, 42.0, 35.0, 14.0, 21.0, 1.0, 20.0, 23.0, 32.0, 45.0, 27.0, 40.0, 27.0, 23.0, 27.0, 20.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]

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
