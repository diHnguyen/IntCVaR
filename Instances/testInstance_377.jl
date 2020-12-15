edge = [1 4; 1 7; 1 11; 1 13; 2 4; 2 5; 2 7; 2 11; 3 2; 3 5; 3 6; 3 7; 3 15; 4 7; 4 11; 4 14; 5 7; 5 9; 5 10; 5 15; 6 8; 6 11; 7 2; 7 6; 7 12; 7 13; 7 14; 8 7; 8 10; 8 11; 8 13; 8 14; 8 15; 9 6; 9 8; 9 13; 9 14; 10 7; 10 15; 11 2; 11 4; 11 13; 11 15; 12 3; 12 4; 12 5; 12 7; 12 11; 13 2; 13 5; 14 7; 14 9; 14 11]
cL_orig = [22.0, 266.0, 19.0, 420.0, 14.0, 23.0, 197.0, 64.0, 3.0, 8.0, 133.0, 151.0, 192.0, 38.0, 82.0, 254.0, 38.0, 148.0, 223.0, 438.0, 37.0, 40.0, 11.0, 13.0, 191.0, 85.0, 74.0, 30.0, 5.0, 88.0, 190.0, 300.0, 48.0, 102.0, 11.0, 132.0, 224.0, 94.0, 145.0, 105.0, 79.0, 21.0, 55.0, 7.0, 10.0, 321.0, 215.0, 13.0, 51.0, 281.0, 112.0, 70.0, 71.0]
cU_orig = [22.0, 266.0, 19.0, 420.0, 14.0, 23.0, 197.0, 94.0, 3.0, 26.0, 133.0, 151.0, 192.0, 102.0, 82.0, 538.0, 38.0, 148.0, 223.0, 438.0, 37.0, 40.0, 17.0, 13.0, 191.0, 85.0, 572.0, 30.0, 5.0, 88.0, 190.0, 300.0, 48.0, 102.0, 11.0, 132.0, 224.0, 94.0, 145.0, 567.0, 79.0, 51.0, 77.0, 7.0, 212.0, 321.0, 215.0, 13.0, 943.0, 281.0, 112.0, 70.0, 71.0]
d = [26.0, 40.0, 20.0, 8.0, 47.0, 3.0, 16.0, 16.0, 24.0, 6.0, 19.0, 14.0, 8.0, 31.0, 16.0, 49.0, 44.0, 15.0, 41.0, 30.0, 45.0, 11.0, 13.0, 14.0, 9.0, 24.0, 17.0, 27.0, 42.0, 11.0, 34.0, 45.0, 7.0, 14.0, 15.0, 7.0, 49.0, 13.0, 36.0, 40.0, 21.0, 15.0, 46.0, 11.0, 27.0, 46.0, 27.0, 5.0, 34.0, 23.0, 28.0, 46.0, 45.0]
Len = length(d)

yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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