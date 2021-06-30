edge = [1 2; 1 11; 1 12; 1 16; 2 3; 2 4; 2 6; 2 9; 2 10; 2 11; 2 13; 2 14; 2 19; 3 2; 3 6; 3 18; 4 2; 4 7; 4 15; 5 4; 5 7; 5 11; 5 16; 5 20; 6 7; 6 8; 6 9; 6 16; 6 19; 7 6; 7 11; 7 15; 7 18; 8 2; 8 7; 8 9; 8 15; 8 19; 9 3; 10 4; 10 6; 10 17; 10 20; 11 6; 11 15; 11 17; 12 4; 12 6; 12 7; 12 16; 12 19; 12 20; 13 4; 13 6; 13 10; 13 12; 13 15; 14 6; 14 7; 14 9; 14 10; 14 18; 14 19; 15 4; 15 5; 15 7; 15 9; 15 13; 16 5; 16 6; 16 8; 16 14; 17 3; 17 8; 17 11; 17 12; 17 20; 18 6; 18 10; 18 19; 19 2; 19 12; 19 13; 19 14; 19 16; 19 18]
cL_orig = [5.0, 82.0, 73.0, 289.0, 8.0, 72.0, 40.0, 101.0, 232.0, 375.0, 30.0, 2.0, 92.0, 3.0, 138.0, 472.0, 35.0, 118.0, 144.0, 48.0, 66.0, 27.0, 26.0, 224.0, 14.0, 73.0, 26.0, 380.0, 214.0, 14.0, 49.0, 113.0, 213.0, 199.0, 8.0, 15.0, 110.0, 327.0, 43.0, 219.0, 129.0, 17.0, 39.0, 41.0, 59.0, 28.0, 127.0, 22.0, 7.0, 95.0, 108.0, 142.0, 352.0, 230.0, 17.0, 11.0, 57.0, 50.0, 43.0, 92.0, 114.0, 47.0, 67.0, 547.0, 103.0, 122.0, 237.0, 11.0, 359.0, 121.0, 297.0, 9.0, 94.0, 191.0, 90.0, 21.0, 78.0, 400.0, 156.0, 5.0, 504.0, 174.0, 263.0, 52.0, 46.0, 21.0]
cU_orig = [15.0, 858.0, 73.0, 289.0, 8.0, 72.0, 110.0, 101.0, 232.0, 375.0, 64.0, 10.0, 1326.0, 3.0, 138.0, 612.0, 59.0, 118.0, 204.0, 48.0, 116.0, 81.0, 26.0, 224.0, 14.0, 73.0, 26.0, 380.0, 214.0, 14.0, 57.0, 113.0, 213.0, 199.0, 8.0, 15.0, 110.0, 327.0, 43.0, 219.0, 129.0, 17.0, 39.0, 411.0, 59.0, 28.0, 265.0, 22.0, 19.0, 95.0, 276.0, 148.0, 352.0, 230.0, 31.0, 11.0, 57.0, 160.0, 303.0, 134.0, 114.0, 47.0, 67.0, 547.0, 103.0, 122.0, 237.0, 67.0, 359.0, 121.0, 297.0, 33.0, 94.0, 191.0, 90.0, 21.0, 78.0, 400.0, 156.0, 5.0, 504.0, 210.0, 263.0, 312.0, 46.0, 69.0]
d = [12.0, 42.0, 37.0, 14.0, 48.0, 19.0, 1.0, 6.0, 48.0, 48.0, 27.0, 28.0, 16.0, 1.0, 8.0, 3.0, 36.0, 39.0, 4.0, 7.0, 36.0, 39.0, 33.0, 31.0, 19.0, 26.0, 9.0, 46.0, 49.0, 31.0, 38.0, 8.0, 25.0, 5.0, 8.0, 12.0, 21.0, 44.0, 38.0, 48.0, 46.0, 40.0, 44.0, 30.0, 38.0, 16.0, 28.0, 8.0, 37.0, 35.0, 4.0, 42.0, 5.0, 6.0, 11.0, 30.0, 2.0, 50.0, 29.0, 3.0, 5.0, 7.0, 19.0, 1.0, 30.0, 38.0, 39.0, 46.0, 32.0, 39.0, 5.0, 49.0, 16.0, 20.0, 19.0, 47.0, 15.0, 38.0, 28.0, 15.0, 49.0, 48.0, 2.0, 47.0, 13.0, 36.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
delta2 = 5
last_node = maximum(edge)
