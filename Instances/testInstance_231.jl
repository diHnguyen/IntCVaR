edge = [1 3; 1 9; 2 3; 2 8; 2 11; 3 2; 3 6; 3 7; 3 13; 4 3; 4 6; 4 11; 4 13; 4 15; 5 3; 5 4; 5 10; 5 14; 5 15; 6 3; 6 8; 6 10; 6 11; 6 14; 7 2; 7 3; 7 8; 7 13; 7 14; 8 2; 8 4; 9 3; 9 6; 9 7; 10 5; 10 9; 10 14; 10 15; 11 5; 11 7; 11 9; 11 10; 11 12; 11 15; 12 3; 12 5; 12 6; 12 14; 13 3; 13 4; 13 7; 13 11; 13 15; 14 2; 14 5; 14 6; 14 7; 14 8; 14 11]
cL_orig = [65.0, 290.0, 19.0, 212.0, 39.0, 33.0, 8.0, 44.0, 45.0, 35.0, 100.0, 27.0, 88.0, 431.0, 28.0, 35.0, 133.0, 253.0, 2.0, 118.0, 90.0, 23.0, 250.0, 2.0, 112.0, 103.0, 20.0, 7.0, 50.0, 42.0, 128.0, 1.0, 3.0, 89.0, 232.0, 41.0, 5.0, 177.0, 290.0, 29.0, 33.0, 19.0, 15.0, 2.0, 16.0, 10.0, 68.0, 34.0, 488.0, 267.0, 37.0, 27.0, 56.0, 32.0, 306.0, 221.0, 45.0, 192.0, 31.0]
cU_orig = [65.0, 290.0, 19.0, 212.0, 39.0, 33.0, 108.0, 44.0, 45.0, 39.0, 100.0, 27.0, 88.0, 431.0, 28.0, 49.0, 133.0, 253.0, 614.0, 118.0, 90.0, 131.0, 250.0, 2.0, 112.0, 103.0, 20.0, 7.0, 50.0, 108.0, 150.0, 1.0, 109.0, 89.0, 232.0, 41.0, 325.0, 305.0, 290.0, 33.0, 77.0, 19.0, 15.0, 64.0, 16.0, 10.0, 360.0, 34.0, 488.0, 267.0, 37.0, 157.0, 56.0, 98.0, 306.0, 221.0, 203.0, 192.0, 31.0]
d = [14.0, 28.0, 50.0, 15.0, 4.0, 50.0, 3.0, 25.0, 31.0, 50.0, 10.0, 33.0, 5.0, 3.0, 17.0, 13.0, 39.0, 46.0, 8.0, 24.0, 7.0, 10.0, 21.0, 35.0, 13.0, 39.0, 9.0, 46.0, 28.0, 35.0, 39.0, 28.0, 9.0, 15.0, 21.0, 2.0, 42.0, 43.0, 17.0, 3.0, 48.0, 12.0, 23.0, 40.0, 9.0, 32.0, 46.0, 23.0, 30.0, 14.0, 5.0, 37.0, 48.0, 7.0, 46.0, 39.0, 6.0, 26.0, 45.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]

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
