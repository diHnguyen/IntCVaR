edge = [1 3; 1 4; 1 5; 1 12; 1 13; 1 14; 1 15; 1 19; 1 20; 2 3; 2 6; 2 9; 2 12; 2 17; 3 10; 3 12; 3 15; 4 9; 4 10; 4 11; 4 15; 5 12; 5 14; 6 2; 6 5; 6 8; 6 11; 6 14; 6 18; 7 4; 8 7; 8 9; 8 10; 8 11; 8 15; 8 16; 9 5; 9 13; 10 3; 10 7; 10 9; 10 16; 10 20; 11 3; 11 8; 11 9; 11 10; 11 12; 11 13; 11 14; 12 6; 12 7; 12 10; 12 14; 12 16; 12 17; 13 2; 13 4; 13 6; 13 7; 13 11; 13 18; 14 9; 14 10; 14 15; 14 20; 15 3; 15 5; 15 6; 15 8; 15 11; 15 20; 16 6; 16 10; 16 11; 16 12; 16 15; 16 19; 17 15; 17 20; 18 6; 18 7; 18 12; 18 20; 19 10; 19 12; 19 15; 1 21; 2 23; 4 25; 6 24; 8 21; 8 23; 13 22; 14 22; 14 24; 15 24; 19 22; 19 23; 19 24; 21 8; 21 10; 21 13; 21 17; 21 18; 21 19; 21 25; 22 2; 22 9; 23 3; 23 6; 23 18; 23 24; 23 25; 24 2; 24 4; 24 7; 24 10; 24 11; 24 14; 24 23; 25 7; 25 11; 25 23]
cL_orig = [24.0, 35.0, 39.0, 109.0, 119.0, 128.0, 138.0, 179.0, 189.0, 6.0, 37.0, 69.0, 88.0, 153.0, 75.0, 88.0, 124.0, 46.0, 49.0, 62.0, 105.0, 65.0, 91.0, 41.0, 3.0, 15.0, 48.0, 74.0, 123.0, 32.0, 10.0, 5.0, 14.0, 33.0, 72.0, 78.0, 41.0, 32.0, 75.0, 35.0, 5.0, 55.0, 89.0, 82.0, 27.0, 22.0, 8.0, 14.0, 12.0, 33.0, 57.0, 52.0, 16.0, 21.0, 39.0, 45.0, 106.0, 81.0, 69.0, 62.0, 23.0, 51.0, 53.0, 28.0, 12.0, 61.0, 120.0, 104.0, 90.0, 71.0, 38.0, 39.0, 99.0, 65.0, 46.0, 37.0, 8.0, 34.0, 25.0, 34.0, 120.0, 112.0, 58.0, 21.0, 95.0, 73.0, 37.0, 200.0, 210.0, 210.0, 183.0, 134.0, 145.0, 92.0, 82.0, 101.0, 94.0, 34.0, 42.0, 49.0, 129.0, 107.0, 83.0, 44.0, 31.0, 22.0, 40.0, 203.0, 125.0, 202.0, 174.0, 45.0, 9.0, 24.0, 215.0, 199.0, 172.0, 138.0, 127.0, 105.0, 10.0, 179.0, 143.0, 19.0]
cU_orig = [24.0, 35.0, 39.0, 109.0, 119.0, 138.0, 138.0, 179.0, 189.0, 6.0, 37.0, 69.0, 106.0, 153.0, 75.0, 88.0, 124.0, 46.0, 65.0, 82.0, 105.0, 65.0, 91.0, 41.0, 23.0, 15.0, 48.0, 84.0, 123.0, 32.0, 10.0, 19.0, 32.0, 33.0, 72.0, 92.0, 41.0, 44.0, 75.0, 35.0, 5.0, 55.0, 109.0, 82.0, 27.0, 22.0, 8.0, 14.0, 30.0, 33.0, 57.0, 52.0, 16.0, 21.0, 39.0, 45.0, 122.0, 91.0, 69.0, 62.0, 23.0, 51.0, 53.0, 42.0, 12.0, 61.0, 120.0, 104.0, 90.0, 71.0, 38.0, 51.0, 99.0, 65.0, 56.0, 37.0, 8.0, 34.0, 25.0, 34.0, 120.0, 112.0, 58.0, 21.0, 95.0, 73.0, 37.0, 200.0, 210.0, 210.0, 183.0, 134.0, 145.0, 92.0, 82.0, 101.0, 94.0, 34.0, 42.0, 49.0, 129.0, 107.0, 83.0, 44.0, 31.0, 22.0, 40.0, 203.0, 125.0, 202.0, 174.0, 45.0, 9.0, 24.0, 215.0, 199.0, 172.0, 138.0, 127.0, 105.0, 10.0, 179.0, 143.0, 19.0]
d = [3.0, 14.0, 8.0, 12.0, 11.0, 5.0, 15.0, 10.0, 3.0, 6.0, 16.0, 16.0, 12.0, 11.0, 12.0, 20.0, 5.0, 16.0, 7.0, 1.0, 3.0, 3.0, 19.0, 14.0, 14.0, 9.0, 8.0, 8.0, 7.0, 16.0, 2.0, 13.0, 5.0, 19.0, 7.0, 11.0, 7.0, 17.0, 11.0, 14.0, 12.0, 17.0, 7.0, 11.0, 1.0, 19.0, 3.0, 18.0, 6.0, 3.0, 1.0, 13.0, 13.0, 7.0, 5.0, 5.0, 3.0, 7.0, 3.0, 17.0, 7.0, 17.0, 18.0, 14.0, 16.0, 3.0, 10.0, 11.0, 2.0, 13.0, 8.0, 2.0, 16.0, 11.0, 6.0, 16.0, 14.0, 8.0, 9.0, 14.0, 18.0, 7.0, 12.0, 19.0, 3.0, 4.0, 11.0, 15.0, 12.0, 14.0, 5.0, 10.0, 15.0, 14.0, 5.0, 10.0, 4.0, 6.0, 15.0, 6.0, 11.0, 16.0, 3.0, 8.0, 12.0, 19.0, 12.0, 9.0, 16.0, 9.0, 6.0, 9.0, 4.0, 6.0, 1.0, 16.0, 18.0, 1.0, 4.0, 18.0, 14.0, 3.0, 10.0, 4.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)
p = [1.0]
g = [SP_init]
h = [0.0]
origin = 1
destination = 20
last_node = maximum(edge)
all_nodes = collect(1:last_node)
M_orig = zeros(Len)
for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end
delta1 = 1.0
delta2 = 5.0
b = 2
last_node = maximum(edge)
Grp = "+5Nodes"
