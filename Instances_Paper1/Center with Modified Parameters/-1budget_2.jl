edge = [1 3; 1 4; 1 15; 1 18; 2 6; 2 14; 2 18; 3 2; 3 5; 3 10; 3 11; 3 14; 3 20; 4 7; 4 8; 4 11; 4 18; 4 19; 5 6; 5 7; 5 9; 5 13; 5 19; 5 20; 6 4; 6 8; 6 14; 7 4; 7 10; 7 12; 7 14; 7 16; 7 20; 8 3; 8 13; 8 19; 8 20; 9 20; 10 4; 10 6; 10 9; 10 15; 10 17; 10 19; 11 6; 11 8; 11 15; 12 2; 12 3; 12 6; 12 11; 12 15; 12 18; 12 19; 12 20; 13 8; 13 15; 14 3; 14 13; 14 19; 14 20; 15 4; 15 5; 15 6; 15 16; 16 4; 16 20; 17 19; 18 5; 18 13; 18 17; 18 20; 19 4; 19 8; 19 10; 19 13; 19 14; 19 15; 19 16; 19 17; 19 20]
cL_orig = [25.0, 32.0, 140.0, 167.0, 36.0, 121.0, 157.0, 7.0, 21.0, 65.0, 73.0, 109.0, 172.0, 22.0, 43.0, 71.0, 138.0, 147.0, 4.0, 23.0, 36.0, 82.0, 140.0, 142.0, 22.0, 22.0, 84.0, 34.0, 29.0, 47.0, 67.0, 89.0, 119.0, 48.0, 50.0, 104.0, 118.0, 109.0, 65.0, 31.0, 5.0, 49.0, 67.0, 95.0, 45.0, 33.0, 38.0, 101.0, 86.0, 58.0, 7.0, 28.0, 49.0, 61.0, 81.0, 47.0, 23.0, 98.0, 3.0, 52.0, 65.0, 106.0, 103.0, 82.0, 8.0, 116.0, 40.0, 19.0, 135.0, 48.0, 15.0, 21.0, 152.0, 109.0, 91.0, 61.0, 49.0, 39.0, 35.0, 15.0, 7.0]
cU_orig = [25.0, 32.0, 140.0, 167.0, 36.0, 121.0, 157.0, 7.0, 21.0, 65.0, 93.0, 119.0, 172.0, 42.0, 43.0, 71.0, 138.0, 147.0, 22.0, 23.0, 36.0, 82.0, 140.0, 156.0, 22.0, 22.0, 84.0, 34.0, 29.0, 57.0, 67.0, 89.0, 137.0, 48.0, 50.0, 114.0, 118.0, 109.0, 65.0, 47.0, 19.0, 49.0, 77.0, 95.0, 45.0, 33.0, 38.0, 101.0, 86.0, 58.0, 7.0, 28.0, 69.0, 71.0, 81.0, 47.0, 23.0, 118.0, 23.0, 52.0, 65.0, 106.0, 103.0, 92.0, 8.0, 116.0, 40.0, 19.0, 135.0, 48.0, 15.0, 21.0, 152.0, 109.0, 91.0, 61.0, 49.0, 39.0, 35.0, 15.0, 7.0]
d = [10.0, 19.0, 11.0, 4.0, 6.0, 7.0, 10.0, 5.0, 18.0, 7.0, 2.0, 20.0, 1.0, 12.0, 15.0, 8.0, 16.0, 14.0, 10.0, 3.0, 14.0, 10.0, 20.0, 4.0, 4.0, 11.0, 14.0, 13.0, 14.0, 20.0, 19.0, 5.0, 10.0, 16.0, 16.0, 11.0, 1.0, 9.0, 3.0, 3.0, 9.0, 8.0, 2.0, 20.0, 6.0, 6.0, 6.0, 11.0, 13.0, 8.0, 7.0, 7.0, 16.0, 19.0, 6.0, 3.0, 12.0, 6.0, 1.0, 10.0, 10.0, 19.0, 7.0, 16.0, 12.0, 1.0, 1.0, 18.0, 7.0, 4.0, 13.0, 12.0, 2.0, 14.0, 3.0, 11.0, 18.0, 16.0, 14.0, 18.0, 1.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
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
delta2 = 4.5
b = 1
last_node = maximum(edge)
Grp = "-1budget"
