edge = [1 3; 1 4; 1 10; 1 12; 1 15; 2 6; 2 9; 2 13; 2 14; 3 5; 3 12; 3 13; 3 14; 4 3; 4 7; 4 13; 5 3; 5 8; 5 9; 5 11; 5 13; 5 14; 6 2; 6 10; 7 3; 7 13; 7 14; 7 15; 8 7; 8 10; 8 11; 8 14; 9 4; 9 5; 9 6; 9 10; 9 12; 9 13; 10 2; 10 4; 10 5; 10 13; 10 14; 10 15; 11 2; 11 4; 11 7; 11 12; 11 13; 12 3; 12 9; 12 13; 13 8; 13 10; 13 11; 14 8; 14 9; 14 12]
cL_orig = [84.0, 116.0, 146.0, 164.0, 686.0, 60.0, 329.0, 112.0, 122.0, 1.0, 58.0, 406.0, 105.0, 8.0, 20.0, 428.0, 16.0, 50.0, 69.0, 25.0, 116.0, 336.0, 8.0, 128.0, 176.0, 169.0, 319.0, 12.0, 17.0, 89.0, 143.0, 253.0, 221.0, 50.0, 46.0, 23.0, 92.0, 65.0, 72.0, 111.0, 24.0, 140.0, 30.0, 145.0, 190.0, 136.0, 134.0, 47.0, 32.0, 322.0, 0.0, 16.0, 36.0, 115.0, 100.0, 10.0, 162.0, 0.0]
cU_orig = [84.0, 116.0, 146.0, 506.0, 686.0, 60.0, 329.0, 780.0, 876.0, 17.0, 644.0, 406.0, 105.0, 70.0, 20.0, 428.0, 16.0, 50.0, 69.0, 25.0, 382.0, 336.0, 52.0, 128.0, 176.0, 169.0, 319.0, 26.0, 17.0, 89.0, 143.0, 253.0, 221.0, 62.0, 186.0, 23.0, 92.0, 207.0, 72.0, 111.0, 24.0, 140.0, 340.0, 145.0, 190.0, 136.0, 134.0, 47.0, 114.0, 322.0, 254.0, 16.0, 36.0, 115.0, 100.0, 176.0, 162.0, 24.0]
d = [2.0, 2.0, 8.0, 15.0, 44.0, 13.0, 50.0, 19.0, 48.0, 18.0, 10.0, 14.0, 26.0, 15.0, 26.0, 26.0, 47.0, 4.0, 46.0, 49.0, 47.0, 46.0, 46.0, 8.0, 21.0, 7.0, 6.0, 24.0, 14.0, 8.0, 42.0, 4.0, 31.0, 4.0, 27.0, 48.0, 9.0, 25.0, 2.0, 33.0, 2.0, 25.0, 6.0, 49.0, 2.0, 41.0, 27.0, 35.0, 1.0, 20.0, 22.0, 25.0, 50.0, 7.0, 39.0, 25.0, 3.0, 42.0]
Len = length(d)

yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
