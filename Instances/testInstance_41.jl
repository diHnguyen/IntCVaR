edge = [1 2; 1 3; 1 4; 1 7; 1 11; 1 12; 2 6; 2 9; 2 11; 2 13; 3 2; 3 13; 3 15; 4 2; 4 3; 4 5; 4 6; 4 7; 5 9; 5 14; 6 2; 6 3; 6 9; 6 11; 6 15; 7 10; 7 11; 8 2; 8 3; 8 5; 8 6; 8 7; 8 10; 8 11; 8 12; 8 15; 9 6; 9 13; 9 15; 10 4; 10 5; 10 11; 11 5; 11 6; 11 9; 11 15; 12 3; 12 5; 12 9; 12 10; 12 11; 12 14; 13 7; 13 9; 13 10; 14 3; 14 4; 14 10; 14 11; 14 13]
cL_orig = [5.0, 22.0, 7.0, 83.0, 453.0, 3.0, 157.0, 251.0, 50.0, 69.0, 16.0, 351.0, 531.0, 94.0, 7.0, 23.0, 12.0, 6.0, 84.0, 203.0, 90.0, 21.0, 115.0, 194.0, 116.0, 44.0, 183.0, 69.0, 11.0, 34.0, 83.0, 19.0, 3.0, 86.0, 39.0, 111.0, 2.0, 87.0, 77.0, 246.0, 173.0, 17.0, 91.0, 76.0, 9.0, 24.0, 88.0, 108.0, 69.0, 8.0, 38.0, 8.0, 219.0, 26.0, 17.0, 27.0, 243.0, 101.0, 86.0, 10.0]
cU_orig = [5.0, 22.0, 7.0, 439.0, 453.0, 3.0, 157.0, 251.0, 390.0, 69.0, 16.0, 351.0, 531.0, 94.0, 7.0, 23.0, 76.0, 98.0, 84.0, 203.0, 90.0, 125.0, 115.0, 194.0, 116.0, 44.0, 183.0, 357.0, 11.0, 34.0, 83.0, 19.0, 3.0, 152.0, 39.0, 373.0, 2.0, 87.0, 317.0, 246.0, 173.0, 17.0, 91.0, 76.0, 9.0, 240.0, 88.0, 108.0, 75.0, 8.0, 38.0, 8.0, 219.0, 124.0, 35.0, 533.0, 243.0, 101.0, 86.0, 64.0]
d = [2.0, 23.0, 32.0, 38.0, 31.0, 5.0, 9.0, 37.0, 1.0, 7.0, 9.0, 24.0, 15.0, 36.0, 30.0, 16.0, 37.0, 47.0, 31.0, 33.0, 17.0, 10.0, 8.0, 4.0, 43.0, 30.0, 43.0, 37.0, 20.0, 7.0, 16.0, 39.0, 27.0, 38.0, 23.0, 34.0, 45.0, 41.0, 27.0, 22.0, 4.0, 30.0, 45.0, 33.0, 14.0, 37.0, 33.0, 44.0, 24.0, 42.0, 12.0, 26.0, 21.0, 40.0, 31.0, 38.0, 21.0, 38.0, 25.0, 2.0]
Len = length(d)

yy = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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