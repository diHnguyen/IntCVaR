edge = [1 2; 1 3; 1 7; 1 8; 1 9; 1 13; 1 14; 2 4; 2 5; 2 6; 2 7; 2 12; 3 6; 3 13; 3 14; 4 3; 4 8; 5 3; 5 7; 5 8; 5 10; 6 5; 6 7; 6 14; 7 2; 7 8; 7 9; 7 12; 8 2; 8 7; 8 10; 9 2; 9 14; 10 4; 10 5; 10 11; 10 12; 10 15; 11 3; 11 5; 11 8; 11 12; 12 11; 12 13; 12 15; 13 2; 13 3; 13 7; 13 9; 13 11; 14 4; 14 5; 14 9; 14 11; 14 12]
cL_orig = [45.0, 25.0, 23.0, 9.0, 182.0, 72.0, 150.0, 11.0, 78.0, 5.0, 49.0, 364.0, 43.0, 47.0, 397.0, 12.0, 14.0, 30.0, 76.0, 13.0, 114.0, 0.0, 18.0, 139.0, 57.0, 20.0, 10.0, 60.0, 88.0, 37.0, 83.0, 312.0, 61.0, 7.0, 44.0, 39.0, 0.0, 82.0, 208.0, 61.0, 26.0, 14.0, 20.0, 28.0, 5.0, 441.0, 20.0, 14.0, 70.0, 56.0, 424.0, 203.0, 226.0, 107.0, 69.0]
cU_orig = [45.0, 25.0, 149.0, 9.0, 182.0, 450.0, 444.0, 153.0, 78.0, 5.0, 145.0, 364.0, 43.0, 227.0, 397.0, 12.0, 14.0, 30.0, 76.0, 13.0, 114.0, 48.0, 18.0, 139.0, 87.0, 30.0, 54.0, 60.0, 428.0, 37.0, 83.0, 312.0, 61.0, 7.0, 44.0, 39.0, 6.0, 82.0, 208.0, 61.0, 26.0, 14.0, 20.0, 28.0, 5.0, 441.0, 20.0, 14.0, 70.0, 56.0, 424.0, 203.0, 226.0, 107.0, 69.0]
d = [18.0, 30.0, 50.0, 12.0, 7.0, 14.0, 11.0, 44.0, 5.0, 14.0, 30.0, 24.0, 43.0, 7.0, 10.0, 27.0, 20.0, 11.0, 49.0, 2.0, 12.0, 14.0, 35.0, 49.0, 49.0, 31.0, 31.0, 9.0, 39.0, 23.0, 36.0, 2.0, 29.0, 13.0, 18.0, 20.0, 28.0, 35.0, 17.0, 35.0, 21.0, 46.0, 10.0, 9.0, 40.0, 2.0, 36.0, 18.0, 20.0, 40.0, 46.0, 45.0, 19.0, 28.0, 8.0]
Len = length(d)

yy = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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