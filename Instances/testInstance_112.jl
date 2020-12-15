edge = [1 3; 1 4; 1 6; 1 12; 2 4; 2 5; 2 14; 3 4; 3 7; 3 13; 3 14; 3 15; 4 3; 4 14; 5 2; 5 4; 5 7; 5 9; 5 10; 5 11; 5 14; 5 15; 6 2; 6 4; 6 9; 6 11; 7 4; 7 5; 7 6; 7 9; 7 10; 7 12; 8 4; 8 6; 8 7; 8 10; 8 13; 9 8; 9 11; 9 12; 9 14; 10 3; 10 11; 10 12; 10 14; 11 6; 11 10; 11 12; 11 13; 11 15; 12 10; 12 11; 12 13; 12 14; 12 15; 13 2; 13 3; 13 5; 14 7; 14 9; 14 15]
cL_orig = [7.0, 23.0, 191.0, 28.0, 13.0, 85.0, 7.0, 12.0, 24.0, 251.0, 160.0, 377.0, 13.0, 37.0, 108.0, 31.0, 48.0, 86.0, 238.0, 114.0, 7.0, 233.0, 7.0, 67.0, 7.0, 208.0, 56.0, 3.0, 32.0, 56.0, 44.0, 9.0, 2.0, 83.0, 16.0, 99.0, 139.0, 13.0, 72.0, 99.0, 29.0, 123.0, 13.0, 77.0, 127.0, 186.0, 9.0, 20.0, 6.0, 33.0, 91.0, 40.0, 16.0, 15.0, 51.0, 211.0, 10.0, 30.0, 45.0, 104.0, 7.0]
cU_orig = [7.0, 23.0, 191.0, 42.0, 35.0, 85.0, 7.0, 12.0, 24.0, 251.0, 208.0, 377.0, 13.0, 37.0, 114.0, 31.0, 48.0, 86.0, 238.0, 114.0, 23.0, 233.0, 163.0, 67.0, 7.0, 208.0, 56.0, 35.0, 32.0, 56.0, 44.0, 95.0, 14.0, 83.0, 16.0, 99.0, 139.0, 13.0, 72.0, 99.0, 29.0, 123.0, 13.0, 77.0, 127.0, 186.0, 9.0, 30.0, 92.0, 97.0, 91.0, 58.0, 16.0, 15.0, 51.0, 613.0, 36.0, 148.0, 591.0, 104.0, 7.0]
d = [21.0, 45.0, 44.0, 11.0, 45.0, 7.0, 29.0, 38.0, 4.0, 35.0, 26.0, 26.0, 8.0, 34.0, 5.0, 21.0, 47.0, 2.0, 11.0, 23.0, 32.0, 14.0, 30.0, 12.0, 21.0, 1.0, 13.0, 36.0, 46.0, 4.0, 33.0, 24.0, 3.0, 35.0, 15.0, 44.0, 17.0, 33.0, 8.0, 10.0, 13.0, 22.0, 14.0, 40.0, 13.0, 47.0, 23.0, 37.0, 8.0, 11.0, 39.0, 44.0, 20.0, 11.0, 5.0, 44.0, 45.0, 1.0, 37.0, 20.0, 19.0]
Len = length(d)

yy = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]

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
