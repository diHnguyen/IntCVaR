edge = [1 5; 1 7; 1 14; 2 4; 2 6; 2 10; 2 11; 3 6; 3 7; 3 9; 3 12; 3 13; 3 14; 4 5; 4 12; 4 13; 4 15; 5 13; 6 3; 6 8; 6 15; 7 8; 7 10; 7 12; 8 2; 8 5; 8 6; 8 13; 9 3; 9 4; 9 11; 9 12; 9 13; 9 14; 10 3; 10 4; 10 8; 10 11; 10 15; 11 2; 11 10; 11 15; 12 4; 12 6; 12 7; 12 8; 12 15; 13 2; 13 6; 13 8; 14 4; 14 5; 14 6; 14 12; 14 13]
cL_orig = [90.0, 136.0, 594.0, 26.0, 103.0, 62.0, 423.0, 133.0, 6.0, 65.0, 194.0, 69.0, 304.0, 49.0, 135.0, 28.0, 480.0, 132.0, 4.0, 86.0, 110.0, 25.0, 34.0, 126.0, 46.0, 146.0, 67.0, 67.0, 283.0, 88.0, 33.0, 24.0, 79.0, 203.0, 36.0, 51.0, 1.0, 46.0, 34.0, 169.0, 0.0, 81.0, 336.0, 193.0, 170.0, 105.0, 132.0, 36.0, 24.0, 104.0, 394.0, 109.0, 214.0, 61.0, 20.0]
cU_orig = [90.0, 136.0, 594.0, 26.0, 103.0, 62.0, 423.0, 133.0, 6.0, 107.0, 194.0, 175.0, 304.0, 49.0, 135.0, 28.0, 480.0, 132.0, 4.0, 86.0, 374.0, 25.0, 34.0, 126.0, 280.0, 146.0, 75.0, 67.0, 283.0, 88.0, 33.0, 24.0, 79.0, 203.0, 36.0, 51.0, 15.0, 46.0, 120.0, 169.0, 20.0, 81.0, 336.0, 193.0, 170.0, 105.0, 132.0, 36.0, 28.0, 104.0, 574.0, 527.0, 562.0, 61.0, 20.0]
d = [28.0, 44.0, 40.0, 43.0, 42.0, 13.0, 37.0, 36.0, 50.0, 44.0, 14.0, 14.0, 1.0, 16.0, 18.0, 48.0, 15.0, 34.0, 8.0, 30.0, 16.0, 24.0, 42.0, 39.0, 19.0, 30.0, 43.0, 41.0, 16.0, 49.0, 13.0, 15.0, 43.0, 15.0, 20.0, 35.0, 13.0, 21.0, 47.0, 30.0, 38.0, 41.0, 28.0, 49.0, 44.0, 27.0, 24.0, 24.0, 2.0, 46.0, 30.0, 14.0, 38.0, 12.0, 32.0]
Len = length(d)

yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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