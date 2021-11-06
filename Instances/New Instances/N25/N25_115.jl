edge = [1 3; 1 4; 1 21; 1 23; 2 10; 2 24; 3 10; 3 15; 3 20; 3 22; 3 25; 4 2; 4 17; 4 19; 4 25; 5 4; 6 5; 6 22; 6 25; 7 3; 7 5; 7 17; 7 20; 7 24; 8 24; 9 5; 9 6; 9 19; 9 22; 10 3; 10 5; 10 7; 11 7; 11 10; 11 14; 11 18; 11 20; 11 24; 12 7; 12 10; 12 11; 12 22; 13 23; 14 2; 14 11; 14 18; 15 14; 15 21; 16 9; 17 11; 18 5; 19 5; 19 6; 19 11; 19 25; 20 15; 21 6; 21 8; 21 9; 21 16; 21 19; 22 16; 22 17; 22 18; 23 6; 23 8; 23 13; 24 7; 24 12]
cL_orig = [22.0, 30.0, 188.0, 223.0, 78.0, 215.0, 71.0, 120.0, 166.0, 187.0, 212.0, 24.0, 128.0, 154.0, 199.0, 13.0, 0.0, 159.0, 185.0, 43.0, 10.0, 95.0, 128.0, 172.0, 153.0, 44.0, 31.0, 92.0, 131.0, 66.0, 51.0, 30.0, 43.0, 10.0, 28.0, 65.0, 91.0, 125.0, 36.0, 17.0, 10.0, 103.0, 100.0, 125.0, 32.0, 40.0, 8.0, 64.0, 70.0, 61.0, 125.0, 144.0, 132.0, 78.0, 59.0, 48.0, 149.0, 128.0, 119.0, 54.0, 18.0, 64.0, 36.0, 38.0, 155.0, 150.0, 95.0, 161.0, 120.0]
cU_orig = [22.0, 30.0, 208.0, 223.0, 78.0, 215.0, 71.0, 120.0, 166.0, 203.0, 232.0, 24.0, 138.0, 154.0, 217.0, 13.0, 14.0, 159.0, 185.0, 43.0, 30.0, 95.0, 128.0, 172.0, 167.0, 44.0, 31.0, 102.0, 131.0, 66.0, 51.0, 30.0, 43.0, 10.0, 28.0, 65.0, 91.0, 139.0, 54.0, 17.0, 10.0, 103.0, 100.0, 125.0, 32.0, 40.0, 8.0, 64.0, 70.0, 61.0, 125.0, 144.0, 132.0, 78.0, 59.0, 48.0, 149.0, 128.0, 119.0, 54.0, 18.0, 64.0, 56.0, 38.0, 175.0, 150.0, 95.0, 173.0, 120.0]
d = [11.0, 10.0, 1.0, 3.0, 6.0, 20.0, 20.0, 5.0, 16.0, 13.0, 17.0, 13.0, 7.0, 10.0, 17.0, 4.0, 18.0, 7.0, 15.0, 4.0, 8.0, 15.0, 11.0, 11.0, 17.0, 5.0, 12.0, 7.0, 3.0, 7.0, 13.0, 20.0, 13.0, 19.0, 11.0, 8.0, 13.0, 13.0, 8.0, 6.0, 8.0, 12.0, 18.0, 17.0, 18.0, 17.0, 4.0, 16.0, 17.0, 11.0, 19.0, 7.0, 20.0, 14.0, 16.0, 5.0, 7.0, 8.0, 16.0, 6.0, 6.0, 7.0, 13.0, 10.0, 8.0, 13.0, 11.0, 16.0, 6.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)
p = [1.0]
g = [SP_init]
h = [0.0]
origin = 1
destination = 25
last_node = maximum(edge)
all_nodes = collect(1:last_node)
M_orig = zeros(Len)
for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end
delta1 = 1.0
delta2 = 2.0
b = 7
last_node = maximum(edge)
β = 0.796