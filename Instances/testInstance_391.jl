edge = [1 5; 1 8; 1 10; 1 13; 1 14; 2 5; 2 8; 2 10; 2 11; 2 12; 2 13; 3 2; 3 4; 3 5; 3 10; 3 13; 3 14; 4 2; 4 3; 4 7; 4 9; 4 13; 4 15; 5 2; 5 3; 5 4; 5 7; 5 11; 5 12; 6 8; 6 9; 6 11; 6 12; 7 4; 7 5; 7 6; 8 2; 8 4; 8 10; 8 11; 9 2; 9 3; 9 4; 9 5; 10 2; 10 5; 10 11; 10 14; 11 2; 11 4; 11 12; 12 2; 12 6; 12 8; 12 13; 12 14; 13 5; 13 6; 13 7; 13 10; 13 11; 13 14; 14 3; 14 13; 14 15]
cL_orig = [15.0, 98.0, 251.0, 566.0, 11.0, 48.0, 174.0, 285.0, 96.0, 73.0, 15.0, 30.0, 27.0, 24.0, 137.0, 479.0, 434.0, 62.0, 5.0, 88.0, 60.0, 60.0, 21.0, 3.0, 88.0, 39.0, 60.0, 28.0, 312.0, 30.0, 75.0, 80.0, 79.0, 125.0, 98.0, 17.0, 28.0, 196.0, 8.0, 16.0, 70.0, 122.0, 44.0, 33.0, 135.0, 121.0, 0.0, 104.0, 272.0, 20.0, 9.0, 26.0, 293.0, 195.0, 23.0, 16.0, 188.0, 310.0, 188.0, 0.0, 0.0, 3.0, 249.0, 39.0, 18.0]
cU_orig = [15.0, 98.0, 551.0, 566.0, 175.0, 48.0, 174.0, 285.0, 96.0, 73.0, 405.0, 30.0, 27.0, 24.0, 137.0, 479.0, 434.0, 62.0, 39.0, 200.0, 60.0, 60.0, 71.0, 3.0, 106.0, 39.0, 60.0, 28.0, 312.0, 30.0, 75.0, 80.0, 237.0, 125.0, 98.0, 17.0, 28.0, 196.0, 152.0, 270.0, 92.0, 122.0, 44.0, 33.0, 135.0, 121.0, 14.0, 104.0, 272.0, 20.0, 9.0, 32.0, 293.0, 195.0, 23.0, 16.0, 456.0, 310.0, 188.0, 286.0, 154.0, 23.0, 249.0, 47.0, 18.0]
d = [6.0, 4.0, 42.0, 24.0, 6.0, 5.0, 35.0, 35.0, 31.0, 34.0, 49.0, 1.0, 7.0, 9.0, 24.0, 35.0, 27.0, 39.0, 48.0, 37.0, 45.0, 49.0, 45.0, 13.0, 49.0, 12.0, 3.0, 9.0, 15.0, 48.0, 32.0, 38.0, 43.0, 41.0, 45.0, 15.0, 35.0, 37.0, 13.0, 6.0, 20.0, 32.0, 16.0, 41.0, 45.0, 3.0, 2.0, 43.0, 30.0, 12.0, 28.0, 35.0, 4.0, 5.0, 33.0, 10.0, 24.0, 7.0, 48.0, 25.0, 26.0, 46.0, 19.0, 30.0, 14.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]

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
