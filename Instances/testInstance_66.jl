edge = [1 3; 1 5; 1 8; 1 9; 1 11; 1 15; 2 3; 2 8; 2 9; 2 10; 2 12; 2 13; 2 14; 3 2; 3 4; 3 7; 3 8; 3 11; 3 14; 3 15; 4 6; 4 7; 4 12; 4 15; 5 3; 5 4; 5 7; 5 9; 5 10; 5 13; 5 15; 6 3; 6 7; 6 8; 6 13; 6 15; 7 2; 7 9; 7 11; 7 12; 7 14; 8 4; 8 5; 8 6; 8 9; 8 10; 8 12; 8 13; 8 14; 8 15; 9 3; 9 8; 9 14; 10 2; 10 5; 10 6; 10 7; 10 8; 10 12; 10 13; 11 3; 11 4; 11 5; 11 6; 11 7; 11 9; 11 10; 11 12; 12 3; 12 5; 12 6; 12 7; 12 10; 12 11; 12 13; 12 15; 13 5; 13 14; 14 5; 14 7]
cL_orig = [27.0, 57.0, 31.0, 56.0, 230.0, 193.0, 9.0, 68.0, 132.0, 56.0, 280.0, 81.0, 165.0, 4.0, 30.0, 82.0, 248.0, 62.0, 92.0, 106.0, 39.0, 147.0, 343.0, 5.0, 2.0, 18.0, 54.0, 31.0, 12.0, 123.0, 396.0, 68.0, 49.0, 36.0, 335.0, 4.0, 33.0, 78.0, 96.0, 129.0, 34.0, 20.0, 125.0, 94.0, 41.0, 85.0, 60.0, 82.0, 165.0, 291.0, 19.0, 27.0, 131.0, 370.0, 22.0, 182.0, 13.0, 42.0, 19.0, 0.0, 114.0, 338.0, 101.0, 127.0, 6.0, 22.0, 10.0, 10.0, 244.0, 54.0, 174.0, 229.0, 75.0, 20.0, 36.0, 90.0, 342.0, 14.0, 372.0, 90.0]
cU_orig = [27.0, 75.0, 31.0, 56.0, 230.0, 193.0, 9.0, 68.0, 132.0, 56.0, 280.0, 81.0, 165.0, 4.0, 56.0, 220.0, 248.0, 62.0, 92.0, 106.0, 103.0, 147.0, 343.0, 107.0, 110.0, 18.0, 54.0, 113.0, 12.0, 643.0, 556.0, 68.0, 49.0, 36.0, 335.0, 4.0, 33.0, 78.0, 96.0, 129.0, 76.0, 86.0, 125.0, 94.0, 41.0, 85.0, 218.0, 266.0, 221.0, 291.0, 341.0, 27.0, 131.0, 370.0, 22.0, 182.0, 251.0, 42.0, 19.0, 6.0, 114.0, 338.0, 101.0, 127.0, 6.0, 46.0, 10.0, 10.0, 628.0, 54.0, 174.0, 229.0, 75.0, 20.0, 62.0, 90.0, 342.0, 14.0, 372.0, 426.0]
d = [10.0, 38.0, 42.0, 31.0, 18.0, 21.0, 37.0, 33.0, 10.0, 35.0, 45.0, 4.0, 24.0, 13.0, 27.0, 34.0, 8.0, 36.0, 22.0, 33.0, 38.0, 4.0, 48.0, 32.0, 26.0, 7.0, 1.0, 19.0, 37.0, 18.0, 6.0, 2.0, 49.0, 14.0, 35.0, 18.0, 21.0, 50.0, 3.0, 14.0, 10.0, 36.0, 34.0, 36.0, 21.0, 22.0, 31.0, 5.0, 14.0, 9.0, 12.0, 12.0, 12.0, 28.0, 28.0, 13.0, 20.0, 6.0, 42.0, 3.0, 1.0, 1.0, 37.0, 26.0, 5.0, 17.0, 21.0, 23.0, 12.0, 30.0, 1.0, 44.0, 15.0, 33.0, 29.0, 17.0, 37.0, 46.0, 33.0, 36.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
