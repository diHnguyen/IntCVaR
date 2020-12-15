edge = [1 2; 1 3; 1 4; 1 7; 2 5; 3 5; 3 7; 3 10; 3 11; 3 13; 3 14; 4 5; 4 7; 4 11; 4 13; 5 4; 5 8; 5 9; 6 3; 6 4; 7 3; 7 4; 7 5; 7 11; 7 12; 7 15; 8 3; 8 5; 8 9; 9 3; 9 5; 9 8; 9 10; 9 11; 10 3; 10 4; 10 5; 10 8; 10 11; 10 15; 11 2; 11 3; 11 4; 12 2; 12 3; 12 4; 13 2; 13 5; 13 6; 13 10; 13 11; 13 12; 13 14; 14 5; 14 8; 14 9; 14 10; 14 11; 14 12]
cL_orig = [8.0, 76.0, 82.0, 34.0, 85.0, 11.0, 1.0, 220.0, 184.0, 214.0, 332.0, 27.0, 101.0, 1.0, 35.0, 2.0, 117.0, 81.0, 84.0, 36.0, 186.0, 19.0, 38.0, 110.0, 35.0, 209.0, 91.0, 147.0, 10.0, 193.0, 55.0, 19.0, 28.0, 95.0, 118.0, 111.0, 135.0, 4.0, 10.0, 66.0, 92.0, 353.0, 77.0, 389.0, 350.0, 77.0, 86.0, 88.0, 85.0, 33.0, 23.0, 30.0, 8.0, 289.0, 58.0, 67.0, 3.0, 78.0, 95.0]
cU_orig = [12.0, 76.0, 90.0, 34.0, 85.0, 11.0, 1.0, 436.0, 184.0, 346.0, 332.0, 27.0, 101.0, 15.0, 35.0, 2.0, 117.0, 81.0, 84.0, 36.0, 186.0, 19.0, 38.0, 188.0, 35.0, 209.0, 91.0, 147.0, 10.0, 193.0, 271.0, 19.0, 28.0, 95.0, 340.0, 111.0, 135.0, 88.0, 10.0, 66.0, 92.0, 353.0, 529.0, 389.0, 528.0, 77.0, 718.0, 88.0, 85.0, 33.0, 131.0, 30.0, 72.0, 289.0, 58.0, 67.0, 5.0, 78.0, 95.0]
d = [16.0, 22.0, 50.0, 43.0, 1.0, 44.0, 13.0, 13.0, 42.0, 17.0, 43.0, 3.0, 35.0, 14.0, 36.0, 10.0, 32.0, 23.0, 23.0, 29.0, 24.0, 31.0, 15.0, 19.0, 36.0, 13.0, 27.0, 1.0, 10.0, 11.0, 15.0, 50.0, 35.0, 46.0, 2.0, 21.0, 33.0, 6.0, 7.0, 36.0, 40.0, 25.0, 28.0, 11.0, 25.0, 21.0, 49.0, 39.0, 4.0, 30.0, 32.0, 32.0, 9.0, 27.0, 24.0, 26.0, 33.0, 48.0, 12.0]
Len = length(d)

yy = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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