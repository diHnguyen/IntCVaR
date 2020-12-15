edge = [1 6; 1 8; 1 9; 1 13; 1 14; 2 3; 2 6; 2 8; 2 9; 2 13; 2 14; 2 15; 3 8; 3 12; 3 14; 3 15; 4 2; 4 6; 4 7; 4 9; 4 10; 4 12; 4 15; 5 3; 5 7; 5 8; 5 10; 5 12; 5 13; 5 15; 6 5; 6 7; 6 11; 7 5; 7 8; 7 9; 7 12; 7 13; 7 15; 8 5; 8 9; 8 11; 8 13; 9 2; 9 6; 9 12; 9 13; 10 2; 10 5; 10 6; 10 9; 10 11; 11 2; 11 3; 11 4; 11 5; 11 7; 11 9; 11 10; 11 12; 11 15; 12 2; 12 4; 12 13; 12 14; 13 4; 13 5; 13 9; 13 14; 14 5; 14 8; 14 9; 14 11]
cL_orig = [190.0, 103.0, 196.0, 54.0, 598.0, 36.0, 70.0, 175.0, 136.0, 247.0, 463.0, 86.0, 112.0, 332.0, 126.0, 381.0, 56.0, 84.0, 58.0, 1.0, 52.0, 381.0, 3.0, 90.0, 26.0, 147.0, 123.0, 145.0, 49.0, 424.0, 10.0, 2.0, 24.0, 25.0, 36.0, 97.0, 25.0, 86.0, 76.0, 92.0, 19.0, 49.0, 30.0, 20.0, 118.0, 80.0, 93.0, 27.0, 31.0, 50.0, 2.0, 46.0, 109.0, 58.0, 96.0, 288.0, 155.0, 74.0, 29.0, 19.0, 75.0, 306.0, 292.0, 6.0, 26.0, 144.0, 138.0, 55.0, 37.0, 347.0, 65.0, 59.0, 128.0]
cU_orig = [190.0, 209.0, 196.0, 812.0, 598.0, 36.0, 70.0, 175.0, 136.0, 247.0, 463.0, 86.0, 112.0, 332.0, 126.0, 381.0, 62.0, 84.0, 166.0, 61.0, 138.0, 381.0, 3.0, 90.0, 32.0, 147.0, 211.0, 223.0, 49.0, 424.0, 38.0, 2.0, 24.0, 25.0, 36.0, 97.0, 25.0, 86.0, 76.0, 92.0, 19.0, 127.0, 156.0, 26.0, 118.0, 80.0, 93.0, 27.0, 191.0, 50.0, 4.0, 46.0, 109.0, 58.0, 96.0, 288.0, 155.0, 74.0, 29.0, 19.0, 75.0, 306.0, 292.0, 6.0, 26.0, 144.0, 138.0, 55.0, 51.0, 511.0, 65.0, 59.0, 128.0]
d = [10.0, 6.0, 32.0, 20.0, 44.0, 36.0, 30.0, 18.0, 26.0, 44.0, 29.0, 7.0, 1.0, 45.0, 21.0, 6.0, 41.0, 14.0, 19.0, 32.0, 39.0, 50.0, 18.0, 7.0, 4.0, 47.0, 18.0, 40.0, 12.0, 43.0, 46.0, 28.0, 34.0, 38.0, 32.0, 43.0, 8.0, 24.0, 26.0, 38.0, 18.0, 14.0, 41.0, 9.0, 13.0, 24.0, 46.0, 38.0, 5.0, 6.0, 36.0, 33.0, 44.0, 16.0, 44.0, 7.0, 12.0, 31.0, 1.0, 45.0, 43.0, 14.0, 48.0, 3.0, 47.0, 34.0, 25.0, 3.0, 8.0, 6.0, 28.0, 44.0, 11.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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