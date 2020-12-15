edge = [1 9; 1 10; 1 12; 1 14; 1 15; 2 5; 2 6; 2 8; 2 11; 3 7; 3 8; 3 12; 3 13; 3 14; 4 3; 4 6; 4 10; 4 14; 5 3; 5 7; 5 8; 5 10; 5 12; 5 13; 5 15; 6 2; 6 3; 6 5; 6 7; 6 12; 7 3; 7 9; 7 13; 7 14; 8 3; 8 12; 9 2; 9 5; 9 8; 9 10; 9 14; 10 5; 10 7; 10 9; 10 11; 10 14; 11 2; 11 3; 11 4; 11 5; 11 6; 11 13; 11 15; 12 4; 12 5; 12 6; 12 15; 13 2; 13 3; 13 4; 13 5; 13 12; 13 15; 14 9; 14 11]
cL_orig = [322.0, 417.0, 41.0, 576.0, 652.0, 131.0, 88.0, 245.0, 49.0, 40.0, 228.0, 423.0, 232.0, 149.0, 19.0, 96.0, 32.0, 493.0, 39.0, 96.0, 91.0, 55.0, 38.0, 81.0, 375.0, 33.0, 141.0, 20.0, 11.0, 11.0, 83.0, 14.0, 3.0, 128.0, 153.0, 127.0, 11.0, 58.0, 2.0, 33.0, 141.0, 107.0, 58.0, 14.0, 10.0, 152.0, 391.0, 185.0, 111.0, 103.0, 219.0, 100.0, 151.0, 64.0, 152.0, 261.0, 65.0, 8.0, 161.0, 148.0, 386.0, 5.0, 75.0, 207.0, 30.0]
cU_orig = [322.0, 417.0, 41.0, 576.0, 652.0, 131.0, 140.0, 245.0, 49.0, 40.0, 228.0, 423.0, 232.0, 319.0, 19.0, 96.0, 460.0, 493.0, 45.0, 96.0, 91.0, 63.0, 38.0, 559.0, 375.0, 329.0, 141.0, 28.0, 71.0, 305.0, 113.0, 14.0, 3.0, 128.0, 153.0, 127.0, 11.0, 58.0, 2.0, 33.0, 141.0, 107.0, 58.0, 14.0, 34.0, 200.0, 391.0, 185.0, 111.0, 131.0, 219.0, 100.0, 151.0, 176.0, 152.0, 261.0, 65.0, 178.0, 161.0, 518.0, 386.0, 27.0, 75.0, 207.0, 30.0]
d = [45.0, 19.0, 25.0, 21.0, 5.0, 27.0, 46.0, 19.0, 2.0, 16.0, 15.0, 46.0, 28.0, 34.0, 7.0, 26.0, 42.0, 27.0, 10.0, 9.0, 50.0, 32.0, 31.0, 8.0, 3.0, 15.0, 15.0, 46.0, 5.0, 24.0, 10.0, 31.0, 37.0, 47.0, 11.0, 22.0, 44.0, 15.0, 36.0, 24.0, 31.0, 15.0, 18.0, 45.0, 34.0, 13.0, 1.0, 37.0, 16.0, 10.0, 8.0, 2.0, 15.0, 38.0, 44.0, 3.0, 43.0, 5.0, 12.0, 6.0, 32.0, 44.0, 48.0, 36.0, 17.0]
Len = length(d)

yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]

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
