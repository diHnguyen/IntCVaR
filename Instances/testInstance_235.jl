edge = [1 9; 1 12; 2 4; 2 5; 2 12; 2 14; 2 15; 3 6; 3 9; 3 10; 3 11; 3 15; 4 3; 4 6; 4 7; 5 7; 5 8; 6 2; 6 7; 6 9; 6 12; 6 15; 7 2; 7 3; 7 5; 7 6; 7 8; 7 10; 7 11; 7 12; 7 14; 8 4; 8 11; 8 12; 8 15; 9 4; 9 6; 9 8; 9 11; 10 5; 10 7; 10 9; 10 11; 10 12; 10 15; 11 6; 11 12; 11 13; 12 4; 12 7; 12 13; 13 3; 13 4; 13 9; 13 10; 14 2; 14 6; 14 8; 14 11; 14 13]
cL_orig = [85.0, 273.0, 17.0, 124.0, 9.0, 111.0, 0.0, 73.0, 9.0, 1.0, 191.0, 358.0, 15.0, 8.0, 22.0, 92.0, 91.0, 5.0, 19.0, 140.0, 257.0, 368.0, 144.0, 18.0, 62.0, 44.0, 46.0, 138.0, 132.0, 40.0, 77.0, 141.0, 133.0, 34.0, 292.0, 15.0, 91.0, 1.0, 24.0, 81.0, 34.0, 0.0, 6.0, 38.0, 205.0, 30.0, 17.0, 25.0, 33.0, 23.0, 37.0, 144.0, 247.0, 135.0, 80.0, 533.0, 27.0, 265.0, 73.0, 11.0]
cU_orig = [85.0, 275.0, 19.0, 124.0, 9.0, 111.0, 278.0, 73.0, 67.0, 1.0, 191.0, 358.0, 15.0, 182.0, 22.0, 92.0, 91.0, 37.0, 19.0, 140.0, 257.0, 368.0, 160.0, 40.0, 62.0, 52.0, 46.0, 138.0, 132.0, 98.0, 77.0, 141.0, 133.0, 236.0, 292.0, 229.0, 183.0, 23.0, 24.0, 291.0, 34.0, 26.0, 6.0, 70.0, 205.0, 30.0, 17.0, 29.0, 33.0, 23.0, 37.0, 154.0, 247.0, 165.0, 80.0, 533.0, 27.0, 265.0, 73.0, 11.0]
d = [29.0, 28.0, 23.0, 47.0, 19.0, 9.0, 25.0, 46.0, 37.0, 14.0, 42.0, 3.0, 49.0, 10.0, 25.0, 10.0, 7.0, 35.0, 47.0, 28.0, 11.0, 41.0, 19.0, 25.0, 3.0, 41.0, 44.0, 27.0, 46.0, 12.0, 10.0, 50.0, 35.0, 1.0, 29.0, 49.0, 49.0, 31.0, 21.0, 25.0, 11.0, 38.0, 5.0, 50.0, 9.0, 46.0, 3.0, 39.0, 37.0, 43.0, 14.0, 35.0, 3.0, 47.0, 46.0, 12.0, 44.0, 35.0, 27.0, 47.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
