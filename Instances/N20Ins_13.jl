edge = [1 2; 1 3; 1 7; 1 8; 1 9; 1 10; 1 12; 1 17; 2 3; 2 5; 2 10; 2 14; 2 15; 2 18; 3 2; 3 6; 3 7; 3 9; 3 17; 3 18; 4 5; 4 8; 4 9; 4 14; 4 20; 5 8; 5 10; 5 19; 5 20; 6 7; 6 9; 6 11; 6 15; 6 16; 6 17; 6 20; 7 5; 7 10; 7 11; 8 3; 8 6; 8 13; 8 16; 8 17; 8 18; 8 19; 9 10; 9 14; 9 19; 10 3; 10 6; 10 7; 10 8; 10 12; 10 13; 10 14; 10 15; 10 16; 10 17; 11 4; 11 5; 11 7; 11 9; 11 13; 11 14; 12 6; 12 8; 12 9; 12 15; 13 3; 13 6; 13 7; 13 10; 13 11; 13 14; 13 19; 14 5; 14 11; 14 16; 14 18; 14 19; 15 2; 15 5; 15 17; 15 18; 16 2; 16 5; 16 6; 16 11; 16 14; 16 17; 16 19; 17 2; 17 14; 17 16; 17 20; 18 5; 18 7; 18 13; 18 15; 19 4; 19 5; 19 10; 19 11; 19 14; 19 15; 19 20]
cL_orig = [48.0, 2.0, 275.0, 68.0, 263.0, 11.0, 254.0, 652.0, 13.0, 5.0, 309.0, 261.0, 64.0, 620.0, 12.0, 36.0, 171.0, 243.0, 670.0, 168.0, 36.0, 12.0, 139.0, 439.0, 423.0, 76.0, 38.0, 93.0, 53.0, 28.0, 38.0, 149.0, 365.0, 493.0, 13.0, 130.0, 32.0, 94.0, 23.0, 150.0, 68.0, 91.0, 234.0, 19.0, 65.0, 185.0, 32.0, 161.0, 128.0, 34.0, 178.0, 2.0, 22.0, 21.0, 63.0, 35.0, 206.0, 34.0, 131.0, 240.0, 214.0, 60.0, 9.0, 57.0, 109.0, 259.0, 164.0, 126.0, 17.0, 68.0, 203.0, 221.0, 67.0, 30.0, 37.0, 133.0, 391.0, 42.0, 82.0, 19.0, 115.0, 602.0, 122.0, 41.0, 1.0, 59.0, 225.0, 19.0, 24.0, 58.0, 47.0, 14.0, 233.0, 113.0, 21.0, 17.0, 407.0, 3.0, 175.0, 143.0, 546.0, 186.0, 391.0, 286.0, 64.0, 44.0, 44.0]
cU_orig = [48.0, 2.0, 275.0, 84.0, 263.0, 11.0, 254.0, 652.0, 13.0, 223.0, 309.0, 261.0, 64.0, 620.0, 12.0, 36.0, 171.0, 243.0, 670.0, 168.0, 36.0, 124.0, 139.0, 439.0, 423.0, 76.0, 38.0, 93.0, 719.0, 28.0, 70.0, 149.0, 365.0, 493.0, 13.0, 130.0, 112.0, 94.0, 57.0, 150.0, 68.0, 91.0, 234.0, 45.0, 325.0, 185.0, 32.0, 161.0, 128.0, 34.0, 178.0, 4.0, 22.0, 21.0, 63.0, 231.0, 206.0, 34.0, 131.0, 240.0, 298.0, 60.0, 9.0, 57.0, 109.0, 259.0, 164.0, 126.0, 17.0, 68.0, 203.0, 221.0, 67.0, 162.0, 37.0, 133.0, 391.0, 42.0, 82.0, 171.0, 115.0, 602.0, 122.0, 41.0, 1.0, 941.0, 225.0, 505.0, 24.0, 58.0, 51.0, 152.0, 233.0, 113.0, 21.0, 17.0, 841.0, 75.0, 175.0, 143.0, 546.0, 364.0, 391.0, 382.0, 64.0, 44.0, 44.0]
d = [6.0, 16.0, 18.0, 38.0, 24.0, 32.0, 44.0, 43.0, 33.0, 41.0, 7.0, 3.0, 16.0, 20.0, 1.0, 49.0, 23.0, 2.0, 11.0, 3.0, 26.0, 2.0, 43.0, 45.0, 16.0, 43.0, 39.0, 20.0, 26.0, 15.0, 15.0, 18.0, 13.0, 29.0, 42.0, 2.0, 4.0, 33.0, 9.0, 46.0, 26.0, 32.0, 9.0, 38.0, 18.0, 32.0, 41.0, 12.0, 32.0, 6.0, 49.0, 18.0, 5.0, 35.0, 34.0, 19.0, 36.0, 19.0, 17.0, 36.0, 46.0, 28.0, 21.0, 31.0, 33.0, 23.0, 25.0, 44.0, 47.0, 38.0, 8.0, 33.0, 9.0, 21.0, 20.0, 45.0, 37.0, 35.0, 42.0, 33.0, 46.0, 10.0, 4.0, 8.0, 49.0, 47.0, 33.0, 40.0, 35.0, 33.0, 10.0, 44.0, 43.0, 42.0, 22.0, 47.0, 9.0, 38.0, 10.0, 41.0, 33.0, 13.0, 34.0, 32.0, 15.0, 24.0, 50.0]
Len = length(d)

yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

c_orig = 0.5*(cL_orig+cU_orig)

SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)


p = [1.0]

g = [SP_init]

h = [0.0]


origin = 1

destination =20

last_node = maximum(edge)
all_nodes = collect(1:last_node)

M_orig = zeros(Len)

for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end

case = 0
delta1 = 1e-6
delta2 = 2
last_node = maximum(edge)
