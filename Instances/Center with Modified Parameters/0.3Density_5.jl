edge = [1 6; 1 10; 1 11; 1 14; 1 16; 1 20; 2 3; 2 5; 2 8; 2 20; 3 10; 4 10; 4 17; 5 2; 5 4; 5 8; 5 13; 5 15; 6 9; 6 17; 7 11; 7 14; 7 18; 8 9; 8 17; 9 3; 9 10; 9 18; 10 3; 10 12; 10 14; 10 18; 11 6; 12 3; 12 4; 12 5; 12 6; 12 10; 12 13; 12 17; 13 2; 13 4; 13 10; 13 12; 14 19; 14 20; 15 2; 15 4; 15 5; 15 7; 15 13; 15 17; 15 20; 16 3; 16 6; 16 7; 16 8; 16 17; 17 15; 17 20; 18 5; 18 6; 18 20; 19 17; 1 5; 1 7; 1 19; 2 7; 2 15; 3 4; 4 12; 5 7; 5 18; 5 20; 6 3; 6 5; 6 13; 6 20; 7 16; 8 6; 8 12; 8 14; 9 2; 9 19; 10 4; 10 9; 10 15; 12 20; 13 7; 13 11; 13 18; 13 19; 14 8; 16 15; 17 2; 17 3; 18 10; 19 4; 19 9; 19 18; 19 20]
cL_orig = [39.0, 89.0, 102.0, 127.0, 155.0, 191.0, 14.0, 31.0, 62.0, 177.0, 63.0, 59.0, 133.0, 32.0, 9.0, 25.0, 80.0, 102.0, 29.0, 114.0, 37.0, 66.0, 102.0, 15.0, 92.0, 55.0, 3.0, 88.0, 74.0, 20.0, 39.0, 75.0, 51.0, 83.0, 73.0, 67.0, 54.0, 22.0, 6.0, 51.0, 114.0, 87.0, 25.0, 7.0, 48.0, 64.0, 129.0, 108.0, 95.0, 78.0, 22.0, 21.0, 48.0, 126.0, 99.0, 87.0, 82.0, 15.0, 25.0, 30.0, 131.0, 115.0, 21.0, 19.0, 40.0, 60.0, 181.0, 52.0, 135.0, 9.0, 79.0, 16.0, 129.0, 148.0, 29.0, 9.0, 73.0, 140.0, 94.0, 18.0, 43.0, 55.0, 69.0, 97.0, 59.0, 8.0, 47.0, 84.0, 59.0, 24.0, 54.0, 57.0, 58.0, 15.0, 145.0, 145.0, 76.0, 155.0, 100.0, 10.0, 5.0]
cU_orig = [53.0, 89.0, 102.0, 141.0, 155.0, 191.0, 14.0, 31.0, 62.0, 177.0, 81.0, 59.0, 133.0, 32.0, 9.0, 25.0, 80.0, 102.0, 39.0, 114.0, 37.0, 82.0, 114.0, 15.0, 92.0, 71.0, 21.0, 102.0, 74.0, 20.0, 51.0, 75.0, 51.0, 93.0, 93.0, 67.0, 64.0, 22.0, 6.0, 51.0, 114.0, 87.0, 25.0, 7.0, 48.0, 64.0, 129.0, 108.0, 95.0, 78.0, 22.0, 21.0, 48.0, 126.0, 99.0, 87.0, 82.0, 15.0, 25.0, 30.0, 131.0, 115.0, 21.0, 19.0, 40.0, 60.0, 181.0, 52.0, 135.0, 9.0, 79.0, 16.0, 129.0, 148.0, 29.0, 9.0, 73.0, 140.0, 94.0, 18.0, 43.0, 55.0, 69.0, 97.0, 59.0, 8.0, 47.0, 84.0, 59.0, 24.0, 54.0, 57.0, 58.0, 15.0, 145.0, 145.0, 76.0, 155.0, 100.0, 10.0, 5.0]
d = [16.0, 2.0, 20.0, 19.0, 1.0, 4.0, 7.0, 19.0, 3.0, 4.0, 12.0, 20.0, 13.0, 11.0, 16.0, 1.0, 17.0, 17.0, 5.0, 17.0, 16.0, 17.0, 3.0, 19.0, 18.0, 10.0, 8.0, 12.0, 14.0, 20.0, 9.0, 11.0, 15.0, 16.0, 10.0, 1.0, 1.0, 6.0, 2.0, 13.0, 18.0, 15.0, 4.0, 8.0, 18.0, 18.0, 8.0, 15.0, 3.0, 11.0, 10.0, 6.0, 4.0, 10.0, 15.0, 17.0, 2.0, 11.0, 7.0, 14.0, 18.0, 11.0, 20.0, 13.0, 12.0, 7.0, 13.0, 8.0, 11.0, 5.0, 5.0, 1.0, 7.0, 16.0, 17.0, 19.0, 18.0, 1.0, 8.0, 19.0, 20.0, 2.0, 16.0, 6.0, 18.0, 8.0, 9.0, 5.0, 5.0, 11.0, 3.0, 17.0, 19.0, 16.0, 4.0, 14.0, 17.0, 18.0, 15.0, 12.0, 13.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
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
Grp = "0.3Density"