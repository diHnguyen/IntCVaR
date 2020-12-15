edge = [1 2; 1 8; 1 12; 1 15; 2 6; 2 8; 3 4; 3 5; 3 6; 3 10; 3 11; 3 15; 4 11; 5 2; 5 6; 5 14; 6 3; 7 6; 7 11; 7 12; 7 13; 8 3; 8 6; 8 7; 8 10; 8 12; 9 2; 9 5; 9 7; 9 12; 9 14; 10 2; 10 4; 10 5; 10 6; 10 7; 10 8; 10 11; 11 2; 11 3; 11 12; 11 13; 11 15; 12 5; 12 8; 12 10; 12 13; 13 4; 13 8; 13 10; 13 11; 13 14; 14 3; 14 10; 14 11; 14 15]
cL_orig = [39.0, 339.0, 68.0, 109.0, 42.0, 54.0, 32.0, 6.0, 88.0, 63.0, 339.0, 213.0, 103.0, 13.0, 30.0, 266.0, 93.0, 24.0, 149.0, 194.0, 204.0, 18.0, 53.0, 4.0, 29.0, 1.0, 279.0, 140.0, 99.0, 112.0, 121.0, 96.0, 102.0, 105.0, 9.0, 10.0, 57.0, 2.0, 20.0, 17.0, 8.0, 75.0, 184.0, 341.0, 19.0, 75.0, 44.0, 198.0, 53.0, 102.0, 27.0, 4.0, 54.0, 98.0, 48.0, 29.0]
cU_orig = [39.0, 339.0, 68.0, 109.0, 42.0, 54.0, 32.0, 12.0, 88.0, 63.0, 339.0, 213.0, 375.0, 87.0, 30.0, 266.0, 93.0, 38.0, 149.0, 194.0, 204.0, 18.0, 53.0, 46.0, 29.0, 9.0, 279.0, 140.0, 99.0, 112.0, 121.0, 504.0, 136.0, 105.0, 9.0, 10.0, 57.0, 2.0, 72.0, 163.0, 8.0, 75.0, 184.0, 341.0, 253.0, 75.0, 44.0, 198.0, 195.0, 102.0, 27.0, 4.0, 510.0, 232.0, 204.0, 29.0]
d = [15.0, 2.0, 33.0, 17.0, 23.0, 41.0, 46.0, 22.0, 34.0, 48.0, 20.0, 34.0, 41.0, 2.0, 32.0, 7.0, 23.0, 36.0, 47.0, 17.0, 4.0, 22.0, 33.0, 39.0, 30.0, 1.0, 12.0, 20.0, 19.0, 14.0, 9.0, 24.0, 20.0, 29.0, 16.0, 40.0, 48.0, 48.0, 30.0, 42.0, 45.0, 15.0, 19.0, 13.0, 40.0, 19.0, 11.0, 36.0, 10.0, 12.0, 46.0, 21.0, 1.0, 36.0, 8.0, 39.0]
Len = length(d)

yy = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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