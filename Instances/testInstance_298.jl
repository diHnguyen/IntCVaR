edge = [1 5; 1 9; 1 10; 1 14; 1 15; 2 11; 2 12; 3 4; 3 6; 3 7; 3 8; 3 10; 3 11; 3 14; 3 15; 4 3; 4 8; 4 9; 4 12; 4 13; 4 15; 5 2; 5 3; 5 4; 5 7; 5 9; 5 12; 6 11; 6 12; 7 2; 7 4; 7 5; 7 6; 7 12; 7 13; 7 15; 8 4; 8 6; 8 7; 8 13; 9 3; 9 5; 9 10; 9 12; 9 14; 10 3; 10 5; 10 8; 10 15; 11 6; 11 10; 11 13; 12 6; 13 3; 13 6; 13 7; 13 9; 13 11; 14 2; 14 9; 14 12; 14 15]
cL_orig = [163.0, 40.0, 426.0, 511.0, 60.0, 347.0, 421.0, 23.0, 128.0, 102.0, 143.0, 112.0, 7.0, 10.0, 383.0, 17.0, 37.0, 107.0, 168.0, 410.0, 119.0, 86.0, 57.0, 23.0, 10.0, 32.0, 20.0, 6.0, 103.0, 52.0, 16.0, 69.0, 39.0, 155.0, 133.0, 107.0, 89.0, 15.0, 32.0, 239.0, 267.0, 39.0, 41.0, 91.0, 82.0, 174.0, 89.0, 21.0, 85.0, 223.0, 6.0, 96.0, 117.0, 148.0, 340.0, 146.0, 26.0, 34.0, 96.0, 18.0, 24.0, 16.0]
cU_orig = [163.0, 40.0, 426.0, 511.0, 60.0, 347.0, 421.0, 43.0, 128.0, 102.0, 143.0, 112.0, 7.0, 936.0, 383.0, 17.0, 37.0, 107.0, 168.0, 410.0, 159.0, 164.0, 57.0, 23.0, 10.0, 32.0, 130.0, 304.0, 103.0, 254.0, 16.0, 69.0, 39.0, 155.0, 133.0, 107.0, 311.0, 73.0, 32.0, 239.0, 267.0, 39.0, 41.0, 171.0, 82.0, 174.0, 89.0, 65.0, 85.0, 223.0, 6.0, 96.0, 117.0, 148.0, 340.0, 250.0, 96.0, 34.0, 128.0, 18.0, 24.0, 16.0]
d = [40.0, 9.0, 2.0, 50.0, 49.0, 14.0, 25.0, 28.0, 28.0, 12.0, 31.0, 28.0, 6.0, 14.0, 49.0, 25.0, 38.0, 41.0, 30.0, 2.0, 18.0, 40.0, 50.0, 21.0, 26.0, 28.0, 36.0, 6.0, 13.0, 4.0, 18.0, 4.0, 22.0, 1.0, 26.0, 18.0, 14.0, 2.0, 19.0, 37.0, 24.0, 1.0, 30.0, 31.0, 28.0, 28.0, 50.0, 23.0, 12.0, 41.0, 47.0, 10.0, 41.0, 41.0, 29.0, 30.0, 37.0, 41.0, 1.0, 29.0, 34.0, 32.0]
Len = length(d)

yy = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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