edge = [1 3; 1 5; 1 6; 1 8; 1 9; 1 15; 2 3; 2 5; 2 6; 2 11; 2 12; 2 13; 2 14; 2 15; 3 2; 3 4; 3 5; 3 6; 3 8; 3 11; 3 12; 3 13; 4 2; 4 5; 4 8; 4 15; 5 2; 5 3; 5 9; 5 10; 5 11; 5 12; 5 13; 6 2; 6 11; 6 12; 6 14; 7 3; 7 4; 7 9; 7 11; 7 12; 7 13; 8 7; 8 10; 8 11; 9 3; 9 4; 9 5; 9 6; 9 7; 9 8; 9 14; 9 15; 10 2; 10 6; 10 11; 10 12; 10 14; 11 14; 11 15; 12 4; 12 5; 12 8; 12 9; 12 14; 13 3; 13 4; 13 7; 13 10; 13 11; 13 15; 14 4; 14 8; 14 11; 14 13]
cL_orig = [53.0, 81.0, 15.0, 28.0, 152.0, 486.0, 17.0, 104.0, 5.0, 66.0, 417.0, 436.0, 554.0, 135.0, 13.0, 5.0, 8.0, 83.0, 21.0, 394.0, 66.0, 368.0, 61.0, 23.0, 23.0, 217.0, 7.0, 28.0, 106.0, 136.0, 211.0, 99.0, 116.0, 179.0, 153.0, 290.0, 248.0, 19.0, 142.0, 43.0, 52.0, 152.0, 125.0, 9.0, 11.0, 134.0, 44.0, 83.0, 109.0, 130.0, 38.0, 38.0, 153.0, 79.0, 259.0, 112.0, 22.0, 7.0, 40.0, 45.0, 193.0, 37.0, 211.0, 177.0, 134.0, 40.0, 168.0, 3.0, 214.0, 33.0, 80.0, 63.0, 457.0, 261.0, 68.0, 41.0]
cU_orig = [53.0, 81.0, 115.0, 28.0, 152.0, 486.0, 35.0, 104.0, 137.0, 66.0, 417.0, 436.0, 554.0, 135.0, 13.0, 5.0, 8.0, 83.0, 21.0, 394.0, 256.0, 368.0, 61.0, 23.0, 29.0, 217.0, 229.0, 34.0, 106.0, 136.0, 211.0, 99.0, 116.0, 179.0, 153.0, 290.0, 430.0, 341.0, 142.0, 43.0, 106.0, 152.0, 125.0, 9.0, 11.0, 134.0, 164.0, 83.0, 109.0, 130.0, 112.0, 38.0, 153.0, 79.0, 259.0, 112.0, 44.0, 7.0, 98.0, 45.0, 193.0, 37.0, 211.0, 177.0, 134.0, 40.0, 168.0, 521.0, 214.0, 141.0, 106.0, 63.0, 457.0, 261.0, 68.0, 41.0]
d = [28.0, 25.0, 43.0, 36.0, 49.0, 9.0, 9.0, 18.0, 44.0, 35.0, 13.0, 17.0, 5.0, 10.0, 10.0, 3.0, 14.0, 40.0, 24.0, 31.0, 2.0, 37.0, 25.0, 3.0, 25.0, 4.0, 28.0, 45.0, 30.0, 30.0, 26.0, 6.0, 43.0, 12.0, 22.0, 5.0, 1.0, 47.0, 22.0, 12.0, 20.0, 40.0, 1.0, 39.0, 11.0, 30.0, 34.0, 37.0, 48.0, 16.0, 6.0, 8.0, 24.0, 43.0, 22.0, 25.0, 39.0, 22.0, 11.0, 17.0, 7.0, 27.0, 16.0, 40.0, 20.0, 17.0, 13.0, 15.0, 28.0, 3.0, 50.0, 28.0, 47.0, 21.0, 2.0, 22.0]
Len = length(d)

yy = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
