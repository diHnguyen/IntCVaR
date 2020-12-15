edge = [1 2; 1 7; 1 9; 1 13; 1 15; 2 3; 2 5; 2 12; 3 4; 3 10; 3 12; 3 14; 4 5; 4 7; 4 9; 4 11; 4 13; 5 2; 5 4; 5 6; 5 7; 5 9; 5 14; 6 7; 6 8; 6 10; 6 13; 6 15; 7 3; 7 9; 7 11; 7 13; 7 15; 8 3; 8 7; 9 2; 9 4; 9 14; 9 15; 10 4; 10 5; 10 6; 10 8; 10 15; 11 2; 11 3; 11 4; 11 6; 11 10; 12 3; 12 5; 12 6; 12 8; 12 9; 12 11; 12 14; 13 2; 13 7; 13 8; 13 9; 13 15; 14 2; 14 4; 14 5; 14 8; 14 12; 14 13]
cL_orig = [50.0, 255.0, 119.0, 155.0, 537.0, 9.0, 130.0, 127.0, 5.0, 22.0, 393.0, 261.0, 12.0, 130.0, 210.0, 78.0, 88.0, 57.0, 27.0, 12.0, 87.0, 144.0, 162.0, 29.0, 22.0, 152.0, 284.0, 422.0, 97.0, 1.0, 102.0, 55.0, 26.0, 23.0, 0.0, 222.0, 152.0, 179.0, 156.0, 256.0, 250.0, 4.0, 66.0, 109.0, 410.0, 309.0, 298.0, 3.0, 3.0, 138.0, 131.0, 29.0, 191.0, 145.0, 12.0, 67.0, 162.0, 56.0, 232.0, 160.0, 60.0, 564.0, 115.0, 102.0, 161.0, 95.0, 12.0]
cU_orig = [50.0, 255.0, 119.0, 155.0, 537.0, 39.0, 130.0, 127.0, 67.0, 22.0, 393.0, 261.0, 12.0, 130.0, 210.0, 126.0, 204.0, 57.0, 27.0, 12.0, 87.0, 144.0, 162.0, 29.0, 22.0, 152.0, 284.0, 422.0, 97.0, 1.0, 206.0, 65.0, 26.0, 23.0, 54.0, 222.0, 198.0, 179.0, 266.0, 256.0, 250.0, 4.0, 66.0, 109.0, 410.0, 309.0, 298.0, 3.0, 3.0, 138.0, 131.0, 29.0, 191.0, 145.0, 12.0, 67.0, 162.0, 56.0, 232.0, 160.0, 60.0, 564.0, 115.0, 102.0, 275.0, 95.0, 12.0]
d = [20.0, 33.0, 46.0, 17.0, 13.0, 6.0, 46.0, 11.0, 32.0, 5.0, 30.0, 30.0, 44.0, 38.0, 21.0, 28.0, 50.0, 25.0, 5.0, 20.0, 1.0, 11.0, 40.0, 7.0, 42.0, 49.0, 13.0, 16.0, 28.0, 37.0, 26.0, 28.0, 35.0, 48.0, 29.0, 47.0, 3.0, 12.0, 37.0, 4.0, 50.0, 4.0, 48.0, 5.0, 21.0, 41.0, 36.0, 19.0, 4.0, 19.0, 19.0, 11.0, 28.0, 20.0, 42.0, 11.0, 11.0, 8.0, 39.0, 9.0, 21.0, 35.0, 41.0, 23.0, 46.0, 24.0, 19.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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