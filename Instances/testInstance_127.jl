edge = [1 6; 1 10; 2 5; 2 6; 2 8; 2 10; 2 13; 2 14; 3 2; 3 4; 3 7; 3 8; 3 9; 3 12; 3 13; 4 2; 4 5; 4 6; 4 8; 4 10; 4 11; 4 13; 4 15; 5 4; 5 7; 5 8; 5 9; 5 11; 5 15; 6 9; 6 11; 6 15; 7 5; 7 9; 7 10; 7 11; 7 13; 8 2; 8 4; 8 7; 8 9; 8 12; 9 2; 9 4; 9 14; 9 15; 10 4; 10 6; 10 7; 10 8; 10 9; 10 12; 10 14; 10 15; 11 3; 11 9; 11 12; 11 13; 12 3; 12 5; 12 6; 12 7; 12 8; 12 9; 13 10; 13 11; 13 14; 14 4; 14 6; 14 7; 14 8; 14 9; 14 11]
cL_orig = [150.0, 185.0, 131.0, 176.0, 72.0, 264.0, 374.0, 175.0, 18.0, 49.0, 147.0, 36.0, 265.0, 429.0, 398.0, 96.0, 37.0, 16.0, 95.0, 54.0, 125.0, 311.0, 506.0, 19.0, 13.0, 50.0, 183.0, 27.0, 380.0, 25.0, 204.0, 27.0, 69.0, 0.0, 48.0, 20.0, 147.0, 235.0, 66.0, 40.0, 42.0, 172.0, 13.0, 169.0, 105.0, 1.0, 289.0, 86.0, 24.0, 11.0, 42.0, 74.0, 2.0, 47.0, 19.0, 98.0, 9.0, 0.0, 167.0, 204.0, 195.0, 183.0, 62.0, 28.0, 102.0, 66.0, 35.0, 361.0, 304.0, 349.0, 129.0, 37.0, 63.0]
cU_orig = [150.0, 185.0, 131.0, 176.0, 72.0, 264.0, 374.0, 175.0, 18.0, 49.0, 147.0, 228.0, 265.0, 429.0, 398.0, 96.0, 37.0, 16.0, 115.0, 54.0, 415.0, 311.0, 506.0, 19.0, 13.0, 50.0, 183.0, 27.0, 380.0, 153.0, 204.0, 639.0, 69.0, 106.0, 48.0, 52.0, 405.0, 235.0, 66.0, 40.0, 42.0, 172.0, 13.0, 169.0, 105.0, 123.0, 289.0, 86.0, 24.0, 11.0, 42.0, 74.0, 2.0, 47.0, 19.0, 98.0, 61.0, 18.0, 167.0, 204.0, 195.0, 237.0, 62.0, 28.0, 102.0, 66.0, 35.0, 361.0, 304.0, 349.0, 129.0, 147.0, 63.0]
d = [49.0, 19.0, 34.0, 44.0, 16.0, 16.0, 12.0, 10.0, 26.0, 2.0, 7.0, 48.0, 45.0, 10.0, 39.0, 32.0, 44.0, 44.0, 20.0, 2.0, 36.0, 30.0, 20.0, 46.0, 47.0, 30.0, 41.0, 23.0, 34.0, 32.0, 11.0, 6.0, 28.0, 8.0, 12.0, 5.0, 42.0, 48.0, 47.0, 10.0, 41.0, 40.0, 10.0, 46.0, 50.0, 22.0, 1.0, 34.0, 49.0, 50.0, 33.0, 40.0, 36.0, 44.0, 12.0, 24.0, 29.0, 12.0, 38.0, 20.0, 24.0, 35.0, 19.0, 31.0, 8.0, 19.0, 46.0, 11.0, 12.0, 6.0, 20.0, 27.0, 36.0]
Len = length(d)

yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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