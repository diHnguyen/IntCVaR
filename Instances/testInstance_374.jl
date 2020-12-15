edge = [1 3; 1 4; 1 6; 1 9; 1 15; 2 3; 2 8; 2 10; 2 12; 2 14; 2 15; 3 2; 3 7; 3 13; 4 3; 4 5; 4 11; 4 15; 5 6; 5 8; 5 9; 5 11; 5 14; 6 2; 6 8; 6 9; 6 13; 6 15; 7 2; 7 6; 7 9; 7 10; 7 14; 8 4; 8 5; 8 6; 8 9; 8 10; 8 13; 9 5; 9 6; 9 7; 9 12; 9 14; 9 15; 10 4; 10 6; 10 12; 11 5; 11 10; 11 12; 12 3; 12 5; 12 10; 12 11; 12 15; 13 2; 13 5; 13 7; 13 9; 13 14; 14 3; 14 5; 14 7; 14 13]
cL_orig = [76.0, 56.0, 112.0, 11.0, 96.0, 13.0, 272.0, 250.0, 298.0, 386.0, 431.0, 50.0, 15.0, 236.0, 39.0, 13.0, 44.0, 240.0, 7.0, 30.0, 63.0, 194.0, 341.0, 187.0, 36.0, 6.0, 276.0, 206.0, 58.0, 7.0, 33.0, 79.0, 132.0, 0.0, 142.0, 94.0, 15.0, 81.0, 28.0, 51.0, 60.0, 68.0, 83.0, 194.0, 141.0, 35.0, 127.0, 2.0, 102.0, 23.0, 11.0, 240.0, 222.0, 45.0, 30.0, 73.0, 399.0, 284.0, 104.0, 52.0, 50.0, 301.0, 84.0, 87.0, 21.0]
cU_orig = [76.0, 56.0, 112.0, 11.0, 96.0, 13.0, 282.0, 250.0, 298.0, 386.0, 431.0, 50.0, 33.0, 236.0, 39.0, 13.0, 322.0, 240.0, 51.0, 30.0, 83.0, 318.0, 341.0, 187.0, 90.0, 6.0, 276.0, 206.0, 58.0, 7.0, 33.0, 79.0, 132.0, 38.0, 142.0, 94.0, 15.0, 81.0, 28.0, 51.0, 60.0, 68.0, 83.0, 194.0, 223.0, 35.0, 127.0, 2.0, 102.0, 23.0, 41.0, 286.0, 222.0, 45.0, 30.0, 73.0, 399.0, 284.0, 104.0, 206.0, 50.0, 429.0, 650.0, 87.0, 21.0]
d = [11.0, 31.0, 22.0, 14.0, 39.0, 45.0, 21.0, 16.0, 23.0, 33.0, 45.0, 43.0, 32.0, 2.0, 12.0, 46.0, 35.0, 11.0, 9.0, 28.0, 11.0, 24.0, 19.0, 47.0, 19.0, 37.0, 14.0, 23.0, 2.0, 8.0, 43.0, 5.0, 48.0, 6.0, 9.0, 21.0, 26.0, 24.0, 48.0, 2.0, 47.0, 46.0, 3.0, 39.0, 42.0, 20.0, 30.0, 20.0, 2.0, 1.0, 48.0, 38.0, 48.0, 50.0, 3.0, 1.0, 6.0, 11.0, 5.0, 25.0, 10.0, 20.0, 45.0, 26.0, 46.0]
Len = length(d)

yy = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
