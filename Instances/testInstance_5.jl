edge = [1 4; 1 5; 1 6; 1 13; 2 3; 2 9; 2 13; 2 14; 3 4; 3 5; 3 8; 3 9; 3 11; 3 12; 3 14; 4 9; 4 15; 5 2; 5 3; 5 4; 5 6; 5 8; 5 12; 5 13; 5 14; 6 2; 6 3; 6 4; 6 5; 6 7; 6 10; 6 13; 6 14; 7 6; 7 9; 7 10; 7 12; 7 13; 8 4; 8 7; 8 10; 8 12; 8 13; 9 10; 10 3; 10 4; 10 8; 10 14; 11 7; 11 14; 11 15; 12 5; 12 6; 12 7; 12 9; 12 11; 12 14; 12 15; 13 2; 13 3; 13 4; 13 5; 13 7; 13 8; 13 11; 14 7; 14 13]
cL_orig = [23.0, 195.0, 70.0, 26.0, 13.0, 143.0, 111.0, 13.0, 26.0, 29.0, 73.0, 25.0, 14.0, 360.0, 383.0, 76.0, 311.0, 143.0, 74.0, 16.0, 43.0, 85.0, 101.0, 226.0, 14.0, 93.0, 39.0, 50.0, 42.0, 23.0, 68.0, 47.0, 112.0, 2.0, 82.0, 140.0, 7.0, 53.0, 168.0, 6.0, 43.0, 3.0, 55.0, 4.0, 345.0, 154.0, 20.0, 134.0, 78.0, 101.0, 17.0, 209.0, 72.0, 47.0, 29.0, 34.0, 19.0, 96.0, 179.0, 240.0, 384.0, 388.0, 178.0, 40.0, 11.0, 242.0, 16.0]
cU_orig = [23.0, 195.0, 70.0, 26.0, 13.0, 143.0, 187.0, 157.0, 36.0, 129.0, 73.0, 25.0, 14.0, 360.0, 597.0, 76.0, 311.0, 143.0, 74.0, 16.0, 43.0, 85.0, 533.0, 226.0, 14.0, 93.0, 39.0, 86.0, 42.0, 23.0, 68.0, 517.0, 112.0, 40.0, 82.0, 140.0, 357.0, 53.0, 168.0, 26.0, 43.0, 3.0, 55.0, 6.0, 345.0, 154.0, 20.0, 134.0, 178.0, 101.0, 263.0, 209.0, 72.0, 275.0, 83.0, 66.0, 19.0, 96.0, 179.0, 240.0, 384.0, 388.0, 228.0, 40.0, 11.0, 338.0, 22.0]
d = [19.0, 17.0, 39.0, 14.0, 38.0, 21.0, 41.0, 5.0, 49.0, 32.0, 26.0, 46.0, 35.0, 9.0, 7.0, 8.0, 30.0, 33.0, 37.0, 41.0, 38.0, 21.0, 39.0, 25.0, 25.0, 31.0, 48.0, 30.0, 43.0, 12.0, 24.0, 50.0, 39.0, 12.0, 38.0, 35.0, 43.0, 9.0, 11.0, 11.0, 8.0, 32.0, 43.0, 49.0, 33.0, 10.0, 13.0, 5.0, 42.0, 41.0, 32.0, 48.0, 17.0, 48.0, 28.0, 22.0, 35.0, 42.0, 45.0, 13.0, 42.0, 46.0, 5.0, 18.0, 45.0, 24.0, 41.0]
Len = length(d)

yy = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0]

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