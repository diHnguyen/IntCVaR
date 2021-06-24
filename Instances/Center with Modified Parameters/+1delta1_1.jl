edge = [1 2; 1 8; 1 9; 1 12; 1 17; 1 20; 2 5; 2 8; 3 4; 3 5; 3 8; 3 12; 3 18; 4 8; 4 11; 4 18; 5 6; 5 7; 5 12; 5 13; 5 15; 5 19; 6 3; 6 19; 7 2; 7 10; 7 14; 7 19; 8 15; 8 18; 9 5; 9 15; 9 20; 10 11; 10 13; 10 14; 11 4; 11 9; 11 18; 12 2; 12 3; 12 4; 12 6; 12 8; 12 17; 12 19; 13 4; 13 8; 13 10; 13 19; 13 20; 14 2; 14 6; 14 7; 14 20; 15 4; 15 8; 15 13; 15 14; 15 19; 16 3; 16 5; 16 8; 16 10; 16 11; 16 14; 16 19; 16 20; 17 2; 17 3; 17 4; 17 15; 17 16; 17 19; 18 13; 18 16; 18 17; 18 20; 19 6; 19 9; 19 15]
cL_orig = [15.0, 66.0, 75.0, 106.0, 157.0, 190.0, 31.0, 60.0, 7.0, 10.0, 48.0, 92.0, 155.0, 44.0, 65.0, 131.0, 7.0, 15.0, 58.0, 82.0, 96.0, 140.0, 28.0, 123.0, 38.0, 27.0, 68.0, 122.0, 69.0, 99.0, 45.0, 55.0, 110.0, 8.0, 25.0, 30.0, 71.0, 11.0, 70.0, 99.0, 91.0, 78.0, 55.0, 43.0, 51.0, 70.0, 86.0, 50.0, 29.0, 59.0, 72.0, 115.0, 73.0, 73.0, 59.0, 106.0, 66.0, 8.0, 11.0, 36.0, 127.0, 105.0, 70.0, 62.0, 36.0, 20.0, 29.0, 39.0, 151.0, 138.0, 126.0, 24.0, 5.0, 12.0, 52.0, 18.0, 13.0, 20.0, 128.0, 104.0, 38.0]
cU_orig = [15.0, 66.0, 75.0, 106.0, 157.0, 190.0, 31.0, 60.0, 7.0, 30.0, 48.0, 92.0, 155.0, 44.0, 65.0, 147.0, 7.0, 15.0, 74.0, 82.0, 96.0, 140.0, 28.0, 135.0, 54.0, 37.0, 68.0, 122.0, 69.0, 99.0, 45.0, 73.0, 110.0, 8.0, 25.0, 46.0, 71.0, 27.0, 70.0, 99.0, 91.0, 78.0, 55.0, 43.0, 51.0, 70.0, 86.0, 50.0, 29.0, 59.0, 72.0, 115.0, 85.0, 73.0, 59.0, 106.0, 66.0, 26.0, 11.0, 36.0, 127.0, 105.0, 88.0, 62.0, 54.0, 20.0, 29.0, 51.0, 151.0, 138.0, 126.0, 24.0, 5.0, 30.0, 52.0, 18.0, 13.0, 20.0, 138.0, 104.0, 38.0]
d = [1.0, 7.0, 20.0, 4.0, 2.0, 10.0, 5.0, 13.0, 4.0, 9.0, 18.0, 1.0, 20.0, 7.0, 11.0, 9.0, 18.0, 12.0, 12.0, 19.0, 6.0, 10.0, 9.0, 15.0, 16.0, 18.0, 9.0, 10.0, 1.0, 16.0, 9.0, 12.0, 1.0, 8.0, 1.0, 3.0, 20.0, 15.0, 19.0, 14.0, 1.0, 9.0, 19.0, 17.0, 7.0, 1.0, 3.0, 9.0, 5.0, 7.0, 4.0, 4.0, 17.0, 6.0, 6.0, 15.0, 18.0, 2.0, 6.0, 14.0, 20.0, 8.0, 20.0, 14.0, 17.0, 19.0, 16.0, 19.0, 14.0, 8.0, 3.0, 16.0, 14.0, 8.0, 2.0, 5.0, 3.0, 2.0, 10.0, 9.0, 12.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
delta1 = 2.0
delta2 = 3.0
b = 2
last_node = maximum(edge)
Grp = "+1delta1"